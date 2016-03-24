#!/usr/bin/env python

"""Run epik on all FDA-approved inhibitors."""

import os
import re
import csv
import traceback

from openmoltools import openeye, schrodinger
from openeye import oechem

def read_molecules(filename):
    """Read a file into an OpenEye molecule (or list of molecules).

    Parameters
    ----------
    filename : str
        The name of the file to read (e.g. mol2, sdf)

    Returns
    -------
    molecule : openeye.oechem.OEMol
        The OEMol molecule read, or a list of molecules if multiple molecules are read.
        If no molecules are read, None is returned.

    """

    ifs = oechem.oemolistream(filename)
    molecules = list()
    for mol in ifs.GetOEMols():
        mol_copy = oechem.OEMol(mol)
        molecules.append(mol_copy)
    ifs.close()

    if len(molecules) == 0:
        return None
    elif len(molecules) == 1:
        return molecules[0]
    else:
        return molecules    

def run_epik(name, filename, residue_name, perceive_bonds=False):
    """Generate conformer with OpenEye omega, protonation states with Schrodinger Epik, and charges with OpenEye AM1-BCC.

    Parameters
    ----------
    name : str
       The name of the output directory to generate.
    filename : str
       The mol2, PDB, or SDF file to read in.
    residue_name : str
       Three uppercase letters to name residue.
    perceive_bonds : bool, optional, default=False
       If True, will use geometry to perceive connectivity.
       This is necessary for PDB files.

    """
    # Generate molecule geometry with OpenEye
    print("Generating molecule %s from %s" % (name, filename))
    oe_molecule = read_molecules(filename)
    if perceive_bonds:
        oechem.OEDetermineConnectivity(oe_molecule)

    # Assign geometry and charges with Omega
    oe_molecule = openeye.get_charges(oe_molecule, max_confs=800, strictStereo=True, normalize=True, keep_confs=1)

    # Create output subfolder
    output_basepath = os.path.join(output_dir, name)
    if not os.path.isdir(output_basepath):
        os.mkdir(output_basepath)
    output_basepath = os.path.join(output_basepath, name)

    # Save mol2 file with residue name = first three uppercase letters
    print "Running epik on molecule {}".format(name)
    mol2_file_path = output_basepath + '-input.mol2'
    residue_name = re.sub('[^A-Za-z]+', '', name.upper())[:3]
    openeye.molecule_to_mol2(oe_molecule, mol2_file_path, residue_name=residue_name)

    # Run epik on mol2 file
    mae_file_path = output_basepath + '-epik.mae'
    schrodinger.run_epik(mol2_file_path, mae_file_path, tautomerize=True,
                         max_structures=32, ph_tolerance=10.0)

    # Convert maestro file to sdf and mol2
    output_sdf_filename = output_basepath + '-epik.sdf'
    output_mol2_filename = output_basepath + '-epik.mol2'
    schrodinger.run_structconvert(mae_file_path, output_sdf_filename)
    schrodinger.run_structconvert(mae_file_path, output_mol2_filename)
    
    # Read SDF file.
    uncharged_molecules = read_molecules(output_sdf_filename)

    # Assign charges.
    charged_molecules = list()
    for (index, molecule) in enumerate(uncharged_molecules):
        print "Charging molecule %d / %d" % (index+1, len(uncharged_molecules))
        try:
            charged_molecule = openeye.get_charges(oe_molecule, max_confs=800, strictStereo=True, normalize=True, keep_confs=None)
            charged_molecules.append(charged_molecule)
        except Exception as e:
            print(e)
            print("Skipping protomer/tautomer because of failed charging.")
        
    # Write molecules.
    charged_mol2_filename = output_basepath + '-epik-charged.mol2'
    ofs = oechem.oemolistream(charged_mol2_filename)
    for (index, molecule) in enumerate(charged_molecules):
        oechem.OEWriteMolecule(ofs, molecule)
    ofs.close()

    # Write state penalites.
    outfile = open('state-penalties.out', 'w')
    for (index, molecule) in enumerate(charged_molecules):
        state_penalty = oechem.OEGetSDDataPairs(molecule, 'r_epik_State_Penalty')
        outfile.write('%16.8f\n' % state_penalty)
    outfile.close()

if __name__ == '__main__':
    input_csv_file = 'clinical-kinase-inhibitors.csv'
    output_dir = 'output'

    # Create output directory
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    # Generate Histidine
    run_epik('HIS', 'ACE-HIS-NME.pdb', 'HS1', perceive_bonds=True)
    run_epik('HIS-HIS', 'ACE-HIS-HIS-NME.pdb', 'HS2', perceive_bonds=True)

