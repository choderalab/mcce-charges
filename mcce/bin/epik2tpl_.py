#!/usr/bin/python
import sys
import re
import math
import datetime
import commands



class ATOM:

    def __init__(self, line):
	#print line
	line = line.split()
	if line[0] != '@<TRIPOS>ATOM':
		self.serial = int(line[0])    # the ID number of the atom at the time the file was created
		self.name =  line[1]           # the name of the atom
		self.xyz = (float(line[2]), float(line[3]), float(line[4]))      #  the x y z coordinate of the atom
		self.SYBYL = line[5]          # the SYBYL atom type for the atom
		self.chainID =  line[6]        # the ID number of the substructure containing the atom
		self.resName = line[7]        # the name of the substructure containing the atom
		self.charge  =  float(line[8]) # the charge associated with the atom
		
		

    '''def setData(self,line):
	line = line.split()
	if line[0] != '@<TRIPOS>ATOM':
		print line[1]
		self.serial = int(line[0])    # the ID number of the atom at the time the file was created
		self.name = line[1]           # the name of the atom
		self.xyz = (float(line[2]), float(line[3]), float(line[4]))      #  the x y z coordinate of the atom
		self.SYBYL = line[5]          # the SYBYL atom type for the atom
		self.chainID = line[6]        # the ID number of the substructure containing the atom
		self.resName = line[7]        # the name of the substructure containing the atom
		self.charge  = float(line[8])'''
		
    def getData(self):
	#return self.__serial
	return self.name
	'''return self.xyz
	return self.SYBYL
	return self.chainID
	return self.resName
	return self.charge'''
		
		
		
    



# This function look for line between DELIMITER1 and DELIMITER2
def GetTheSentences(file, str1, str2):
    start_rx = re.compile(str1)
    end_rx = re.compile(str2)

    start = False
    output = []
    with open(file, 'rb') as datafile:
         for line in datafile.readlines():
             if re.match(start_rx, line):
                 start = True
             elif re.match(end_rx, line):
                 start = False
             if start:
                  output.append(line)
    return output

def loadMol(fname):
    # read into atom lines into atom records
	atoms = GetTheSentences(fname, '@<TRIPOS>ATOM', '@<TRIPOS>BOND')
	atoms_list = []
	for atom in atoms:
		atom_obj = ATOM(atom)
		#print atom_obj.getData()
		atoms_list.append(atom_obj)
	
	
	return atoms_list

def write_natom(tpl,atoms_list):
	print "sgve"
	#tpl.write(atoms_list[0])
	

def write_comment_header(tpl):
    tpl.write('####################################\n')
    tpl.write('# Topology File for:\n')
    tpl.write('# {}\n'.format(tpl))
    tpl.write('#\n')
    tpl.write('# Created on: {}\n'.format(datetime.date.today()))
    tpl.write('#\n')
    tpl.write('# Created with: epik2tpl.py\n')
    tpl.write('####################################\n')
    tpl.write('\n')
    return

	
	
def write_tpl(tpl,fname):
    #atoms_list = loadMol(fname)
	#print atoms_list[0].name
	# Write the MCCE2 .tpl file header.
    write_comment_header(tpl)

	# Write the number of atoms associated with each conformer.
    write_natom(tpl,atoms_list)
	
	



def main(fname):
    '''import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', action="store", help = 'file to parse')
    parser.add_argument('ligand', action="store",
                        help = '3 letter code of ligand to extract from file')
    parser.add_argument('chain_identifier', action='store',
                        help = 'chain id for ligand')
    parser.add_argument('pH', action='store',type=float,
                        help = 'pH for which parameters are to be created.')
    parser.add_argument('-p', action='store_false',dest='ideal',default =True,
                        help = 'Extract ligand from PDB HETATM entries.')
    parser.add_argument('-r', action='store_false',dest='reverse_order', default=False,
                        help = 'print conformers starting from negative ions')
    options = parser.parse_args()'''
	#original = "EXAMPLE"
    removeExt = fname.replace(".mol2", "")
    with open(removeExt+'.tpl','w') as tpl:
       	write_tpl(tpl,fname)
	
if __name__ == "__main__":
    if len(sys.argv) < 2:
       	print "Specify a .mol file"
       	sys.exit()
    #loadMol(sys.argv[1])
    main(sys.argv[1])
	
    
