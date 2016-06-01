#!/usr/bin/python
'''
epik2tpl.py is created by Salah Salah with the help of the tpl_maker.py created by Denise M Kilburg for the Gunner Lab at CUNY-City College for the
generation of ligand topology files for use in MCCE from epik output file from John Chodera
lab.  
epik2tpl.py requires a epik file output .mol2 format.


'''

import sys
import re
import math
import datetime
import commands
import os
import re
import tempfile


class Atom(object):
    def __init__(self,name,idnum,element):
        self.name = name.strip()
        self.coords = []
        self.idnum = idnum
        self.connects = []
        self.element = element.strip()
        self.hybrid = ""

class ATOM():
    serial = ''
    name = ''
    xyz = (0.0, 0.0, 0.0)
    SYBYL = ''
    chainID = ''
    resName = ''
    charge = 0.0
	#hybrid = ''
	
    def __init__(self,line):
		line = line.split()
		if line[0] != '@<TRIPOS>ATOM':
			self.serial = int(line[0])    # the ID number of the atom at the time the file was created
			self.name =  line[1]           # the name of the atom
			self.xyz = (float(line[2]), float(line[3]), float(line[4]))      #  the x y z coordinate of the atom
			self.SYBYL = line[5]          # the SYBYL atom type for the atom (hybrid)
			self.chainID =  line[6]        # the ID number of the substructure containing the atom
			self.resName = line[7]        # the name of the substructure containing the atom
			self.charge  =  float(line[8]) # the charge associated with the atom
			#self.hybrid  = line[5] 
			
    def writeout(self):
        print "%d %4s %3s %c %8.3f%8.3f%8.3f" % \
              (self.serial, self.name, self.resName, self.chainID,
               self.xyz[0], self.xyz[1], self.xyz[2])
			
class Conformer():
    def __init__(self,name,charge, pKa=None, state_penalty=None):
		self.name = name
		self.charge = charge
		self.pKa = pKa
		self.state_penalty = state_penalty		

class Protomers():

	def __init__(self,conformers_charges,conformers_names):
		True_False = [ conformers_charges[i]==conformers_charges[i-1] for i in range(len(conformers_charges)) ]
		True_False.append(True) # True_False is 1 less then conformers_charges
		self.ligand_list = []
		self.conformers_charges_list = []
		if len(conformers_charges) != len(conformers_names):
			a = len(conformers_names)/len(conformers_charges)
			ligand_org = conformers_names[1][0:3]
			#tpl.write('CONFLIST {}        {}BK '.format(ligand_org,ligand_org))
			#self.ligand_list.append(ligand_org)
			number_of_atoms = 0
			counter = 1
			counter_two = 1
			counter_three = 1
			for conf_charge in conformers_charges:
				conformer = Conformer(conformers_names[counter][0:3],conf_charge)
				if int(conformer.charge) > 0:
					conf_charge_ = str(conformer.charge)
				elif int(conformer.charge) < 0:
					conf_charge_ = "-" + str(conformer.charge)
				else:
					conf_charge_ = str(conformer.charge)
				counter = counter + a
				ligand = conformer.name + conf_charge_ + str(counter_two)
				self.ligand_list.append(ligand)
				self.conformers_charges_list.append(conf_charge_)
				if True_False[counter_three] is False:
					counter_two = 0
				counter_two += 1
				counter_three += 1
		else:
			ligand_org = conformers_names[0][0:3]
			#tpl.write('CONFLIST {}        {}BK '.format(ligand_org,ligand_org))
			#self.ligand_list.append(ligand_org)
			counter = 1
			for conf_name in conformers_charges:
				conformer = Conformer(ligand_org,conf_name)#conformers_charges[counter-1])
				if int(conformer.charge) > 0:
					conf_charge_ = str(conformer.charge)
				elif int(conformer.charge) < 0:
					conf_charge_ = "-" + str(conformer.charge)
				else:
					conf_charge_ = str(conformer.charge)
				ligand = conformer.name + conf_charge_ + str(counter)
				self.ligand_list.append(ligand)
				self.conformers_charges_list.append(conf_charge_)
				if True_False[counter] is False:
					counter = 0
				counter += 1
		
		# Determine number of atoms
		
		return 
			

			
			
# This function look for lines between DELIMITER1 and DELIMITER2 and returns the lines + number of sections
def GetTheSentences(file, str1, str2):
    start_rx = re.compile(str1)
    end_rx = re.compile(str2)
    number_of_sections = 0  # This should be number of conformers
    start = False
    output = []
    with open(file, 'rb') as datafile:
         for line in datafile.readlines():
             if re.match(start_rx, line):
                 start = True
		 number_of_sections += 1
             elif re.match(end_rx, line):
                 start = False
             if start:
                  output.append(line)
    return output

def write_conformers(tpl,conformers_charges,conformers_names):
	protonation_states = Protomers(conformers_charges,conformers_names)
	tpl.write('CONFLIST {}        {}BK '.format(protonation_states.ligand_list[0][:3],protonation_states.ligand_list[0][:3]))
	for conf in protonation_states.ligand_list:
		tpl.write('{} '.format(conf))
	tpl.write('{} '.format(protonation_states.ligand_list[0][:3] + "DM"))
	tpl.write('\n')
	tpl.write('\n')
	return 
	
def write_natom(tpl,conformers_charges,conformers_names,numberofatoms):
	protonation_states = Protomers(conformers_charges,conformers_names)
	
	template = '{0:9}{1:11}{2:1}\n'
	tpl.write(template.format('NATOM',protonation_states.ligand_list[0][:3]+"BK", '0'))
	i = 0
	for conf in protonation_states.ligand_list:
		tpl.write(template.format('NATOM',
										conf, numberofatoms[i]))
		i += 1
	tpl.write(template.format('NATOM',protonation_states.ligand_list[0][:3]+"DM", '0'))
	tpl.write('\n')
	return

def write_iatom(tpl,conformers_charges,conformers_names,numberofatoms,atoms_list):
	protonation_states = Protomers(conformers_charges,conformers_names)
	
	template = '{0:9s}{1:6s}{2:>5s}{3:5d}\n'
	count = 0
	count_two = 0
	count_three = 0
	count_four = 0
	
	for conf_name in protonation_states.ligand_list:
		count_four += numberofatoms[count_two] + 1
		#print 'From ' + str(count_three) + ' to ' + str(count_four)
		for atom_name in atoms_list[count_three:count_four]: 
			if len(atom_name.resName) < 1:
				continue
			else:
				tpl.write(template.format('IATOM', conf_name,
									atom_name.name, count))
			count += 1
			count_three += 1
		count_three += 1
		#count_four += 1
		count_two += 1
		count = 0
		tpl.write('\n')
	
def write_ATOMNAME(tpl,conformers_charges,conformers_names,numberofatoms,atoms_list):
	
	protonation_states = Protomers(conformers_charges,conformers_names)
	
	template = '{0:9s}{1:6s}{2:5d}{3:>4s}\n'
	count = 0
	count_two = 0
	count_three = 0
	count_four = 0
	
	for conf_name in protonation_states.ligand_list:
		count_four += numberofatoms[count_two] + 1
		#print 'From ' + str(count_three) + ' to ' + str(count_four)
		for atom_name in atoms_list[count_three:count_four]: 
			if len(atom_name.resName) < 1:
				continue
			else:
				tpl.write(template.format('ATOMNAME', conf_name, \
                                  count,atom_name.name))
			count += 1
			count_three += 1
		count_three += 1
		#count_four += 1
		count_two += 1
		count = 0
		tpl.write('\n')

def write_basic_info(tpl):
	tpl.write('#1.Basic Conformer Information: name, pka, em, rxn.\n')
	tpl.write('#23456789A123456789B123456789C\n')
	return
		
def write_proton(tpl,conformers_charges,conformers_names):

	protonation_states = Protomers(conformers_charges,conformers_names)
	tpl.write('\n# PROTON SECTION: PROTON means charge:\n\n')
	template = '{0:9}{1:11}{2:5}\n'
	i = 0
	for conf in protonation_states.ligand_list:
		proton_comment = '# This should be ' + str(protonation_states.conformers_charges_list[i]) + ' , but for now we will set it to zero\n'
		tpl.write(proton_comment)
		
		tpl.write(template.format("PROTON",conf, \
                                      '0'))
    
		i += 1
	tpl.write(template.format("PROTON",protonation_states.ligand_list[0][:3]+"DM", '0'))
	tpl.write('\n')
	
def write_pka(tpl,conformers_charges,conformers_names,epik_State_Penalty):
	
	protonation_states = Protomers(conformers_charges,conformers_names)
	
	tpl.write('# Solution pKa Section: pKa data from CRC Handbook of Chemistry and Physics\n')
	tpl.write('# pka is set to zero\n')
	template = '{0:9s}{1:6s}     {2:8.3f}\n'
	pka = 0 # We may need to fix this 
	kT = 298 * 6.022e23 * 1.381e-23 / 4184 # kcal/mol for 298 K
	pH = 7.4
	
	i = 0
	for conf in protonation_states.ligand_list:
		conformer = Conformer(conf,protonation_states.conformers_charges_list[i])
		if conformer.pKa == None:
			tpl.write(template.format("PKA", conf, 0.0))
		else:
			tpl.write(template.format("PKA", conf, 0.0))
		i += 1
	tpl.write(template.format("PKA", protonation_states.ligand_list[0][:3]+"DM", 0.0))
	tpl.write('\n')
	
	return
		
def write_electron(tpl,conformers_charges,conformers_names):
	
	protonation_states = Protomers(conformers_charges,conformers_names)
	
	tpl.write("#ELECTRON SECTION:\n")
	template = '{0:9}{1:11}{2:5}\n'
 
	for conf in protonation_states.ligand_list:
		tpl.write(template.format("ELECTRON",conf,"0.0"))
	tpl.write(template.format("ELECTRON",protonation_states.ligand_list[0][:3]+"DM","0.0"))
	tpl.write('\n')
	return

def write_EM(tpl,conformers_charges,conformers_names): 
	
	protonation_states = Protomers(conformers_charges,conformers_names)
	template = '{0:9}{1:11}{2:5}\n'
	tpl.write("# EM SECTION:\n")
	
	for conf in protonation_states.ligand_list:
		tpl.write(template.format("EM",conf,"0.0"))
	tpl.write(template.format("EM",protonation_states.ligand_list[0][:3]+"DM","0.0"))
	tpl.write('\n')
	return

def write_RXN(tpl,conformers_charges,conformers_names): 
	protonation_states = Protomers(conformers_charges,conformers_names)
	template = '{0:9}{1:11}{2:5}\n'
	tpl.write("# REACTION FIELD ENERGY SECTION:\n")
	
	for conf in protonation_states.ligand_list:
		tpl.write(template.format("RXN",conf,"0.0"))
	tpl.write('\n')
	return
	
def write_atom_param_section(tpl,protonation_states,atoms_list,vdw_dict):
	
	a = len(atoms_list)/len(protonation_states.ligand_list)
	tpl.write("# Atom Parameters:\n")
	tpl.write("# Van Der Waals Radii. See source for reference\n")
	template = '{0:9}{1:7}{2:5} {3:7}\n'
	for atom in atoms_list[1:a]:
		element = ''.join([i for i in atom.name if not i.isdigit()])[:1]
		#element = "".join(set(element))
		#element = oechem.OEGetAtomicSymbol(atom.GetAtomicNum()).upper()
		tpl.write(template.format("RADIUS",protonation_states.ligand_list[0][:3], \
                                  atom.name,vdw_dict[element.upper()]))
	
	tpl.write('\n')
	return
	
def write_charges(tpl,protonation_states,atoms_list,numberofatoms):
	
	tpl.write("# Atoms Charges:\n")
	template = '{0:9}{1:7}{2:3} {3:7}\n'
	count = 0
	count_two = 0
	count_three = 0
	count_four = 0
	
	for conf in protonation_states.ligand_list:
		count_four += numberofatoms[count_two] + 1
		for atom_name in atoms_list[count_three:count_four]: 
			if len(atom_name.resName) < 1:
				continue
			else:
				#print atom_name.charge
				tpl.write(template.format("CHARGE",conf, \
                                  atom_name.name, atom_name.charge))
	
				count += 1
				count_three += 1
		count_three += 1
		count_two += 1
		count = 0
		tpl.write('\n')
	
	return
	
def write_bond_angel(tpl,protonation_states,atoms_list,numberofatoms):
	template = '{0:9}{1:11}{2:5}\n'
	
	for conf in protonation_states.ligand_list:
		
		tpl.write(template.format("BOND_ANG",protonation_states.ligand_list[0][0:3],"0.0"))
	tpl.write('\n')
	return

def write_extra(tpl, protonation_states, epik_State_Penalty):
	
	template = '{0:9s}{1:6s}{2:5s}{3:8.3f}\n'
	tpl.write("# EXTRA energy for tautomers:\n")
	i = 0
	for state_penalty in epik_State_Penalty:
		tpl.write(template.format("EXTRA", protonation_states.ligand_list[i], "", float(state_penalty)))
		i += 1
	tpl.write('\n')
	return
	
# This function returns a list of strings with the total charge per @<TRIPOS>ATOM section in the mol2 output
def get_conformers_charges(atoms_list):
	conformers = []
	numberofatoms = []
	# Looking for the start of the conformer (atoms number starts with 1)
	conf_start = []
	count = 0
	for atom in atoms_list:
		if atom.serial == 1:
			conf_start.append(count)
		count += 1
	conf_start = [x-1 for x in conf_start[1:]]
	conf_start = [1] + conf_start
	conf_start.append(len(atoms_list))
	
	
	# Add charges within the same section
	number_of_atoms = 0
	total_charge = 0
	for i in range(0,len(conf_start)-1):
		for j in range(conf_start[i],conf_start[i+1]):
			charge = atoms_list[j].charge
			total_charge += charge
			number_of_atoms+= 1
		#conformers.append("{0:.0f}".format(total_charge))
		
		conformers.append(total_charge)
		if i == 0:
			numberofatoms.append(number_of_atoms)
		else:
			numberofatoms.append(number_of_atoms-1)
		total_charge = 0
		number_of_atoms = 0
		
	
	# Convert conformer charges to int
	int_conformers = [] # This is a list of strings
	for conf in conformers:
		conf_ = "{0:.0f}".format(conf)
		if conf_ == '-0':
			conf_ = '0'
		int_conformers.append(conf_)
		
	return int_conformers, numberofatoms

def get_conformers_names(SUBSTRUCTURE_list):
	conformer_name_list = []
	
	for line in SUBSTRUCTURE_list:
		line = line.split()
		if line[0] != '@<TRIPOS>SUBSTRUCTURE':
			#print "ae"
			conformer_name_list.append(line[1])
		#print line
	#conformer_name_list = sorted(set(conformer_name_list))
	#print conformer_name_list
	return conformer_name_list	

def write_comment_header(tpl):
    tpl.write('####################################\n')
    tpl.write('# Topology File for:\n')
    tpl.write('# {}\n'.format(tpl.name))
    tpl.write('#\n')
    tpl.write('# Created on: {}\n'.format(datetime.date.today()))
    tpl.write('#\n')
    tpl.write('# Created with: epik2tpl.py\n')
    tpl.write('####################################\n')
    tpl.write('\n')
    return

def loadMol(fname):
	atoms = GetTheSentences(fname, '@<TRIPOS>ATOM', '@<TRIPOS>BOND')
	atoms_list = []
	for atom in atoms:
		atoms_obj = ATOM(atom)
		atoms_list.append(atoms_obj)
	
	return atoms_list



bond_table = {'H':'s',
			  'H.1':'s',
			  'C.3':'sp3',
			  'C.2':'sp2',
			  'C.1':'sp',
			  'C.ar':'sp2',
			  'C.cat':'unknown',
			  'N.3':'sp3',
			  'N.2':'sp2',
			  'N.1':'sp',
			  'N.ar':'sp2',
			  'N.4':'sp3',
			  'O.3':'sp3',
			  'O.2':'sp2',
			  'O.co2':'unknown',
			  'O.spc':'unknown',
			  'O.t3p':'unknown',
			  'S.3':'sp3',
			  'S.2':'sp2',
			  'S.O':'sp2',
			  'S.O2':'unknown',
			  'P.3':'sp3',	
			  }

bond_table_2 = {'1':'sp3',
				'2':'sp2',
				'3':'sp',
				'ar':'sp2',
				}			  
			  
bond_table_3 = [
    [ "C", "C", '1', 'sp3'], ## C-C
    [ "C", "C", '2', 'sp2'], ## C=C
    [ "C", "C", '3', 'sp'],  ## C-tb-C, tb = triple bond
    [ "C", "N", '1', 'sp3'], ## C-N
    [ "C", "N", '2', 'sp2'], ## C=N
    [ "C", "N", '3', 'sp' ], ## C-tb-N
    [ "C", "O", '1', 'sp3'],  ## C-O
    [ "C", "O", '2', 'sp2'], ## C=O
    [ "C", "S", '1', 'sp3'], ## C-S
    [ "C", "S", '2', 'sp2'], ## C=S
    [ "C", "S", '3', 'sp'],  ## C-tb-S
    [ "C", "P", '1', 'sp3'],
    [ "C", "BR",'1','sp3'],
    [ "C", "CL",'1', 'sp3'],
	[ "C", "Cl",'1', 'sp3'],
    [ "C", "F", '1', 'sp3'],
    [ "C", "I", '1', 'sp3'],
    [ "C", "X", 'ar', 'sp2'], 
	
	[ "O", "C", '1', 'sp3'],
    [ "O", "C", '2', 'sp2'],
    [ "O", "S", '1', 'sp3'],
    [ "O", "S", '2', 'sp2'],
    [ "O", "P", '1', 'sp3'],
    [ "O", "P", '2', 'sp2'],
    [ "O", "N", '1', 'sp3'],
    [ "O", "N", '2', 'sp2'],
    [ "O", "O", '1', 'sp3'],
	[ "O", "O", '2', 'sp2'],
	[ "O", "O", '3', 'sp'],
    [ "O", "X", 'ar', 'sp'], 
	
	[ "N", "C", '1', 'sp3'],
    [ "N", "C", '2', 'sp2'],
    [ "N", "C", '3', 'sp'],
    [ "N", "N", '1', 'sp3'],
    [ "N", "N", '2', 'sp2'],
    [ "N", "O", '1', 'sp3'],
    [ "N", "O", '2', 'sp2'],
    [ "N", "P", '1', 'sp3'],
    [ "N", "P", '2', 'sp2'],
    [ "N", "S", '1', 'sp3'],
    [ "N", "S", '2', 'sp2'],
    [ "N", "SI",'1', 'sp3'],
	[ "N", "Si",'1', 'sp3'],
	[ "N", "X", 'ar', 'sp2'], 
    
	[ "P", "P", '1', 'sp3'],
    [ "P", "P", '2', 'sp2'],
    [ "P", "N", '1', 'sp3'],
    [ "P", "C", '1', 'sp3'],
    [ "P", "O", '1', 'sp3'],
    [ "P", "S", '1', 'sp3'],
    
	[ "S", "N", '3', 'sp'],
    [ "S", "N", '2', 'sp2'],
    [ "S", "O", '3', 'sp2'],
    [ "S", "O", '2', 'sp'],
    [ "S", "C", '2', 'sp2']
]


			  
def write_con_header(tpl): 
    string = "ires" + " " + "conn "
    tpl.write("#ONNECT" + "   " + "conf" + " " + "atom" + "  " + "orbital" \
                  + "  " +string*4 + "\n")
    tpl.write("#ONNECT" +" " + "|-----|----|---------|" + "----|"*10 + "\n")
    return

def search_bond_table(connectionType,second_string):
	if connectionType in bond_table:
		return bond_table[connectionType]
	elif second_string in bond_table_2:
		return bond_table_2[second_string]
	else:
		return 'unknownH'

def write_con_section(tpl,protonation_states,atoms_list,bond_list,numberofatoms):
	template1 = '{0:9}{1:5} {2:^4} {3:^10}'
	template2 = '{0:^4} {1:^4} '
	a = 0
	b = 0
	my_list = []
	for item in atoms_list:
		if len(item.name) < 1:
			a += 1
			my_list.append(b)
		else:
			b += 1
	my_list.append(b)
	
	counter = 0
	a = 0
	b = len(atoms_list)
	k = 0
	for item in bond_list:
		# atomstolookin = atoms_list[0:66]
		# atomstolookin = atoms_list[67:133]
		# atomstolookin = atoms_list[134:199]
		# atomstolookin = atoms_list[200:265]
		if item == '@<TRIPOS>BOND\n':
			counter += 1
			tpl.write('\n')
			tpl.write('\n')
			# Write header of connection
			tpl.write("#  " + protonation_states.ligand_list[k] + "\n")
			write_con_header(tpl)
			k += 1
			print '=================='
			#continue
		else:
			a += 1
			counter += 1
			fields = item.split()
			#print k
			atomA = str(atoms_list[int(fields[1])+my_list[k-1]+k-1].name)
			second_string = fields[3]
			atomA_hybrid = search_bond_table(str(atoms_list[int(fields[1])+my_list[k-1]+k-1].SYBYL),second_string) #bond_table[str(atoms_list[int(fields[1])+my_list[k-1]+k-1].SYBYL)]
			atomB = str(atoms_list[int(fields[2])+my_list[k-1]+k-1].name)
			
			
			tpl.write(template1.format("CONNECT",protonation_states.ligand_list[k-1],atomA,atomA_hybrid))
			tpl.write(template2.format("0", atomB))
			#tpl.write(template2.format("0", atomB))
			tpl.write('\n')
			
			#print atomA + '[' + atomA_hybrid + ']' + ' connected to ' + atomB
			#print connectionType
	
		
		
		
	
	return	
	


## vdw_dict contains vdw radius of atoms as written in
## Bondi,A. van Der Waals Volumes and Radii,
## J.Phys.Chem, 63, 3, 1964 and http://periodic.lanl.gov/index.shtml
## (Los Alamos National Lab Periodic Table)

vdw_dict = {'H':1.20,
            'He':1.40,
            'C':1.70,
            'N':1.55,
			'B':1.92,  # Added by Salah 
            'O':1.52,
            'F':1.47,
            'NE':1.54,
            'SI':2.10,
            'P':1.80,
            'S':1.80,
            'CL':1.75,
			'Cl':1.75,
            'AR':1.88,
            'AS':1.85,
            'SE':1.90,
            'BR':1.85,
            'KR':2.02,
            'TE':2.06,
            'I':1.98,
            'XE':2.16,
            'ZN':1.39,
            'CU':1.40,
            'HG':1.55,
            'CD':1.58,
            'NI':1.63,
            'PD':1.63,
            'AU':1.66,
            'AG':1.72,
            'MG':1.73,
            'PT':1.75,
            'LI':1.82,
            'U':1.86,
            'GA':1.87,
            'PB':2.02,
            'SN':2.17,
            'NA':2.27,
            'FE':0.00,
            'K':2.75}



#################  write_tpl  ###########################################
def write_tpl(tpl,fname,epik_State_Penalty_file):
	import argparse
	
	#Save epik_State_Penalty into a list
	with open(epik_State_Penalty_file, 'r') as f:
		epik_State_Penalty = [line.strip() for line in f]

    # atoms_list is a list of objects of type ATOM
	atoms_list = loadMol(fname)
    
    # Write the MCCE2 .tpl file header.
	write_comment_header(tpl)
	
	conformers_charges, numberofatoms = get_conformers_charges(atoms_list)
	SUBSTRUCTURE_list = GetTheSentences(fname, '@<TRIPOS>SUBSTRUCTURE', '@<TRIPOS>MOLECULE')
	conformers_names = get_conformers_names(SUBSTRUCTURE_list)
	protonation_states = Protomers(conformers_charges,conformers_names)
	
	
	
	# Generate list of conformers with different protonation and tautomer states.
	write_conformers(tpl,conformers_charges,conformers_names)
	
	# Write block describing number of atoms present in each conformer.
	write_natom(tpl,conformers_charges,conformers_names,numberofatoms)

	# Write block describing name of atoms present in each conformer.
	write_iatom(tpl,conformers_charges,conformers_names,numberofatoms,atoms_list)
	
	# Write block describing name of atoms present in each conformer.
	write_ATOMNAME(tpl,conformers_charges,conformers_names,numberofatoms,atoms_list)
	
	# Write Basic Conformer Information: name, pka, em, rxn.
	write_basic_info(tpl)
	
	# PROTON SECTION: PROTON means charge
	write_proton(tpl,conformers_charges,conformers_names)
	
	# Write pka section  
	write_pka(tpl,conformers_charges,conformers_names,epik_State_Penalty)
	
	# Write ELECTRON section 
	write_electron(tpl,conformers_charges,conformers_names)
	
	# Write EM section 
	write_EM(tpl,conformers_charges,conformers_names)
	
	# Write RXN 
	write_RXN(tpl,conformers_charges,conformers_names)
	
	# connect =======
	'''parser = argparse.ArgumentParser()
	parser.add_argument('whatever', action="store", help = 'file to parse')
	parser.add_argument('-p', action='store_false',dest='ideal',default =True, help = 'specifies file is Protein.pdb')
	parser.add_argument('-r', action='store_false',dest='reverse_order', default=True, help = 'print conformers starting from negative ions')
	options = parser.parse_args()
	
	pdb_file_name = fname.split('-')
	options.ligand = protonation_states.ligand_list[0][:3]
	options.charge = 1
	options.chain_identifier = 'A'''
	
	

	
	
	# con section
	
	bond_list = GetTheSentences(fname, '@<TRIPOS>BOND', '@<TRIPOS>SUBSTRUCTURE')
	
	write_con_section(tpl,protonation_states,atoms_list,bond_list,numberofatoms)
		
	
	
	
	# Atom Parameters
	write_atom_param_section(tpl,protonation_states,atoms_list,vdw_dict)
	
	# Atom Charges
	write_charges(tpl,protonation_states,atoms_list,numberofatoms)
	
	# BOND_ANG
	write_bond_angel(tpl,protonation_states,atoms_list,numberofatoms)
	
	# State Penalty
	write_extra(tpl,protonation_states,epik_State_Penalty)
	
	

	return


#################  Main  ###########################################
def main(fname_base,fname, epik_State_Penalty_file): 
	lines = [line .strip() for line in open(fname).readlines()]
	removeExt = fname_base.replace(".mol2", "")
	three_characters = removeExt[:3]
	paramDir = "param"
	if not os.path.exists(paramDir):
		os.makedirs(paramDir)
	#print 
	if len(lines[1]) > 3:
		lines[1] = three_characters.upper()
		
	with open(os.path.join(paramDir, lines[1]+'.tpl'),'w') as tpl:
		write_tpl(tpl,fname,epik_State_Penalty_file)
	
	

if __name__ == "__main__":
	if len(sys.argv) < 2:
		print "Specify an input\nFor example python epik2tpl.py Bosutinib "
		sys.exit()
	fname_base = sys.argv[1]
	
	
	fname_merged = fname_base+'-merged.mol2'
	fname_state_penalties =  fname_base+'-state-penalties.out'
	
	main(fname_base,fname_merged,fname_state_penalties)

	
	