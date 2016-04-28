#!/usr/bin/env python

"""Run epik on all FDA-approved inhibitors."""

import os
import re
import csv
import traceback
import numpy as np

from openeye import oechem
from openmoltools import openeye, schrodinger

MAX_ENERGY_PENALTY = 10.0 # kT


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

def DumpSDData(mol):
    print ("SD data of", mol.GetTitle())
    #loop over SD data
    for dp in oechem.OEGetSDDataPairs(mol):
        print (dp.GetTag(), ':', dp.GetValue())
    print ()

def enumerate_conformations(name, smiles):
    """Generate geometry and run epik."""
    # Generate molecule geometry with OpenEye
    print "Generating molecule {}".format(name)
    oe_molecule = openeye.smiles_to_oemol(smiles)
    try:
        oe_molecule = openeye.get_charges(oe_molecule, keep_confs=1)
    except RuntimeError as e:
        traceback.print_exc()
        print "Skipping molecule " + name
        return

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
    schrodinger.run_epik(mol2_file_path, mae_file_path, tautomerize=False,
                         max_structures=100, min_probability=np.exp(-MAX_ENERGY_PENALTY), ph=7.4)

    # Convert maestro file to sdf and mol2
    output_sdf_filename = output_basepath + '-epik.sdf'
    output_mol2_filename = output_basepath + '-epik.mol2'
    schrodinger.run_structconvert(mae_file_path, output_sdf_filename)
    schrodinger.run_structconvert(mae_file_path, output_mol2_filename)

    # Read SDF file.
    ifs_sdf = oechem.oemolistream()
    ifs_sdf.SetFormat(oechem.OEFormat_SDF)
    ifs_sdf.open(output_sdf_filename)
    sdf_molecule = oechem.OEMol()
    uncharged_molecules = read_molecules(output_sdf_filename)

    # Read MOL2 file.
    ifs_mol2 = oechem.oemolistream()
    ifs_mol2.open(output_mol2_filename)
    mol2_molecule = oechem.OEMol()
    uncharged_molecules = read_molecules(output_sdf_filename)

    # Assign charges.
    charged_molecules = list()
    index = 0
    while oechem.OEReadMolecule(ifs_sdf, sdf_molecule):
        molecule = oechem.OEReadMolecule(ifs_mol2, mol2_molecule)
        index += 1
        print "Charging molecule %d / %d" % (index, len(uncharged_molecules))
        try:
            # Charge molecule.
            charged_molecule = openeye.get_charges(sdf_molecule, max_confs=800, strictStereo=False, normalize=True, keep_confs=None)

            # Store tags.
            oechem.OECopySDData(charged_molecule, sdf_molecule)

            charged_molecules.append(charged_molecule)
        except Exception as e:
            print(e)
            print("Skipping protomer/tautomer because of failed charging.")

    # Clean up
    ifs_sdf.close()
    ifs_mol2.close()

    # Write molecules.
    charged_mol2_filename = output_basepath + '-epik-charged.mol2'
    ofs = oechem.oemolostream(charged_mol2_filename)
    for (index, charged_molecule) in enumerate(charged_molecules):
        oechem.OEWriteMolecule(ofs, charged_molecule)
    ofs.close()

    # Write first molecule as PDB
    charged_pdb_filename = output_basepath + '-epik-charged.pdb'
    ofs = oechem.oemolostream(charged_pdb_filename)
    for (index, charged_molecule) in enumerate(charged_molecules):
        if index > 0: break
        oechem.OEWriteMolecule(ofs, charged_molecule)
    ofs.close()

    # Write state penalites.
    outfile = open(output_basepath + '-state-penalties.out', 'w')
    for (index, charged_molecule) in enumerate(charged_molecules):

        # Get Epik data.
        epik_Ionization_Penalty = float(oechem.OEGetSDData(charged_molecule, "r_epik_Ionization_Penalty"))
        epik_Ionization_Penalty_Charging = float(oechem.OEGetSDData(charged_molecule, "r_epik_Ionization_Penalty_Charging"))
        epik_Ionization_Penalty_Neutral = float(oechem.OEGetSDData(charged_molecule, "r_epik_Ionization_Penalty_Neutral"))
        epik_State_Penalty = float(oechem.OEGetSDData(charged_molecule, "r_epik_State_Penalty"))
        epik_Tot_Q = int(oechem.OEGetSDData(charged_molecule, "i_epik_Tot_Q"))

        outfile.write('%16.8f\n' % epik_State_Penalty)
    outfile.close()

if __name__ == '__main__':
    input_csv_file = 'clinical-kinase-inhibitors.csv'
    output_dir = 'output'

    # Create output directory
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    # Parse csv file
    with open(input_csv_file, 'r') as csv_file:
        for name, smiles in csv.reader(csv_file):
            enumerate_conformations(name, smiles)

    # Generate Histidine
    #enumerate_conformations('Histidine', 'O=C([C@H](CC1=CNC=N1)N)O')
