#!/usr/bin/env python

"""Run epik on all FDA-approved inhibitors."""

import os
import re
import csv
import traceback
import numpy as np

from openmoltools import openeye, schrodinger

MAX_ENERGY_PENALTY = 10.0 # kT

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
    schrodinger.run_structconvert(mae_file_path, output_basepath + '-epik.sdf')
    schrodinger.run_structconvert(mae_file_path, output_basepath + '-epik.mol2')


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
    enumerate_conformations('Histidine', 'O=C([C@H](CC1=CNC=N1)N)O')
