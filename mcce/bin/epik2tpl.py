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
	
    def __init__(self,line):
		line = line.split()
		if line[0] != '@<TRIPOS>ATOM':
			self.serial = int(line[0])    # the ID number of the atom at the time the file was created
			self.name =  line[1]           # the name of the atom
			self.xyz = (float(line[2]), float(line[3]), float(line[4]))      #  the x y z coordinate of the atom
			self.SYBYL = line[5]          # the SYBYL atom type for the atom
			self.chainID =  line[6]        # the ID number of the substructure containing the atom
			self.resName = line[7]        # the name of the substructure containing the atom
			self.charge  =  float(line[8]) # the charge associated with the atom
			
			
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
		print True_False
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
					conf_charge_ = "+" + str(counter_two)
				elif int(conformer.charge) < 0:
					conf_charge_ = "-" + str(counter_two)
				else:
					conf_charge_ = str(conformer.charge)
				counter = counter + a
				ligand = conformer.name + conf_charge_ #+ str(counter_two)
				print ligand
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
				conformer = Conformer(ligand_org,conf_name)
				if int(conformer.charge) > 0:
					#conf_charge_ = "+" + str(conformer.charge)
					conf_charge_ = "+" + str(counter)
				elif int(conformer.charge) < 0:
					#conf_charge_ = "-" + str(conformer.charge)
					conf_charge_ = "-" + str(counter)
				else:
					conf_charge_ = "0" + str(counter) #str(conformer.charge)
				
				ligand = conformer.name + conf_charge_ #+ str(counter)
				print ligand
				self.ligand_list.append(ligand)
				self.conformers_charges_list.append(conf_charge_)
				print counter
				if True_False[counter] is False:
					counter = 0
				counter += 1
		
		# Determine number of atoms
		
		return 
			
class Pdb(object):
    def add_atom(self,atom):
        if atom.idnum > self.max_idnum:
            self.max_idnum = atom.idnum
        self.atombyid[atom.idnum] = atom
        self.atom_list.append(atom)

    def add_hetatm(self,line,options):
        if options.ligand == line[17:20] \
                and options.chain_identifier == line[21]:
            atom = Atom(line[12:16],int(line[6:11].strip()),line[76:78])
            atom.coords.append(float(line[30:38]))
            atom.coords.append(float(line[38:46]))
            atom.coords.append(float(line[46:54]))
            self.add_atom(atom)

    def add_hydrogen(self,connect):
        self.max_idnum += 1
        atom = Atom("H"+str(self.hydrogen_count),self.max_idnum,'H')
        self.hydrogen_count += 1
        atom.connects.append(connect)
        atom.hybrid = 's'
        self.add_atom(atom)
        return atom

    def add_conect(self,line):
        idnum = int(line[6:11])
        if self.atombyid.has_key(idnum):
            atom = self.atombyid[idnum]
            for i in [11,16,21,26]:
                c = line[i:i + 5].strip()
                if c != '':
                    connect_id = int(c)
                    if self.atombyid.has_key(connect_id):
                        atom.connects.append(connect_id)
        return
        
    def read_pdb(self,options,file):
        for line in file:
            if line.startswith("CONECT"):
                self.add_conect(line)
            elif line.startswith("HETATM"):
                self.add_hetatm(line,options)
        
    def read_ideal_pdb(self,options,file):
        for line in file:
            if line.startswith("ATOM  "):
                self.add_hetatm(line,options)
            elif line.startswith("CONECT"):
                self.add_conect(line)

    def read(self,options):
        with open(options.filename, 'r') as file:
            if options.ideal:
                self.read_ideal_pdb(options,file)
            else:
                self.read_pdb(options,file)

    def __init__(self,options):
        self.atom_list = []
        self.atombyid = {}
        self.max_idnum = 0
        self.hydrogen_count = 1
        self.read(options)
			
			
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
	

	
'''def write_con_section(tpl,protonation_states,atoms_list):

	i = 0
	for conf in protonation_states.ligand_list:
		tpl.write("#  " + conf + "\n")
		write_con_header(tpl)
		template1 = '{0:9}{1:5} {2:^4} {3:^10}'
		template2 = '{0:^4} {1:^4} '
		hydrogen_skips = []
		j = protonation_states.conformers_charges_list[i]
		for atom in atoms_list:
			hyb_table = ['s', 'sp', 'sp2', 'sp3', 'sp3d', 'sp3d2']
			#hybrid = hyb_table[oechem.OEGetHybridization(atom)]
		
	return'''

	
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

def mk_conformers(options,conf_chrage):
    conformers = []
    charge = options.charge
    #count = options.charge
    while charge < 0:
        con = options.ligand+str(charge)
        con = Conformer(con,charge)
        conformers.append(con)
        charge += 1
    con = conf_chrage
    con = Conformer(con,charge)
    conformers.append(con)
    return conformers

def mk_vectors(atom1,atom2,atom3):
    vect1 = []
    vect2 = []
    for i in range(3):
        vect1.append(float(atom2.coords[i])-float(atom1.coords[i]))
        vect2.append(float(atom3.coords[i])-float(atom1.coords[i]))
    return vect1,vect2

def mk_dot_prod(vect1,vect2):
    dot_prod = 0
    for i in range(3):
        dot_prod = dot_prod + vect1[i]*vect2[i]
    return dot_prod

def magnitude(vector):
    mag = math.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)
    return mag

def find_angle(vect1,vect2):
    try:
        cos_ang = mk_dot_prod(vect1,vect2)/(magnitude(vect1)*magnitude(vect2))
        angle = math.degrees(math.acos(cos_ang))
        return angle
    except:
        return None

def write_con_header(tpl): 
    string = "ires" + " " + "conn "
    tpl.write("#ONNECT" + "   " + "conf" + " " + "atom" + "  " + "orbital" \
                  + "  " +string*4 + "\n")
    tpl.write("#ONNECT" +" " + "|-----|----|---------|" + "----|"*10 + "\n")
    return
		
def write_con_section(options,tpl,pdb,conformers):
    conf_start = 0
    conf_end = len(conformers) + 1
    #print conf_end
    for conformer in conformers:
        tpl.write("#  " + conformer.name + "\n")
        write_con_header(tpl)
        template1 = '{0:9}{1:5} {2:^4} {3:^10}'
        template2 = '{0:^4} {1:^4} '
        hydrogen_skips = []
        j = conformer.charge
        #print "====="
        for atom in pdb.atom_list:
            #print atom.name
            if not(atom.idnum in hydrogen_skips):
                tpl.write(template1.format("CONNECT",conformer.name, \
                                           atom.name,atom.hybrid))
                if j < 0 and atom.element == 'O' and atom.hybrid == 'sp3':
                    connects = atom.connects[:-1]
                    hydrogen_skips.append(connects[-1])
                    j += 1
                else:
                    connects = atom.connects
                for cid in connects:
                    tpl.write(template2.format("0", pdb.atombyid[cid].name))
                tpl.write('\n')
        tpl.write('\n')
    tpl.write('\n')
    return

def o_connected(atom,pdb):
    for idnum in atom.connects:
        if pdb.atombyid[idnum].element == 'O':
            return True
    return False
	
def search_bond_table(a1,a2,distance):
    for row in bond_table:
        if row[0] == a1.element and row[1] == a2.element \
                and row[2] <= distance and distance <= row[3]:
            return row[4]
    return 'unknownB'
	
def search_hyd2_table(a1,connects):
    for row in hydrogen_table2:
        if row[0] == a1 and row[1] == connects:
            return row[2]
    return 'unknownH'
	
def add_hybrids(pdb):
    for atom in pdb.atom_list:
        connects = len(atom.connects)
        #print atom.element + '---' +str(connects)
        atom.hybrid = search_hyd2_table(atom.element,connects)
    return
	
def search_hydrogen_table(atom):
    for row in hydrogen_table:
        if row[0] == atom.element and row[2] == atom.hybrid \
                and row[1] == len(atom.connects):
            return row[3]
    return 0  	
	
def add_hydrogens(pdb):
    for atom1 in pdb.atom_list:
        con_len = len(atom1.connects)
        if con_len >= 2:
            atom2 = pdb.atombyid[atom1.connects[0]]
            atom3 = pdb.atombyid[atom1.connects[1]]
            (vector1, vector2) = mk_vectors(atom1,atom2,atom3)
            angle = find_angle(vector1,vector2)
            atom1.hybrid = geometry(angle)
        elif con_len == 1 and atom1.element != 'H':
            atom2 = pdb.atombyid[atom1.connects[0]]
            bond_len = bond_length(atom1,atom2)
            orbital = search_bond_table(atom1,atom2,bond_len)
            atom1.hybrid = orbital

        num_hydrogens = search_hydrogen_table(atom1)
        for i in range(num_hydrogens):
            h_atom = pdb.add_hydrogen(atom1.idnum)
            atom1.connects.append(h_atom.idnum)
	
def bond_length(atom1,atom2):
    vector = []
    for i in range(3):
        vector.append(atom1.coords[i] - atom2.coords[i])
    length = magnitude(vector)
    return length	

# Returns hybridization based on angles. 
def geometry(angle):
    if angle > 106 and angle < 115:
        return 'sp3'
    if angle < 183 and angle > 177:
        return 'sp'
    if angle < 125 and angle > 115:
        return 'sp2'
    else:
        return 'unknownA'
	
#################  Main  ###########################################

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


## [Atom, number of atoms connected, hybridization, number of hydrogens to add]
hydrogen_table = [
    ['C', 4, 'sp3', 0],
    ['C', 3, 'sp3', 1],
    ['C', 3, 'sp2', 0],
    ['C', 2, 'sp3', 2],
    ['C', 2, 'sp2', 1],
    ['C', 2, 'sp', 0],
    ['C', 1, 'sp3', 3],
    ['C', 1, 'sp2', 2],
    ['C', 1, 'sp', 1],
    ['N', 3, 'sp3', 0],
    ['N', 3, 'sp2', 0],
    ['N', 2, 'sp3', 1],
    ['N', 2, 'sp2', 1],
    ['N', 2, 'sp', 0],
    ['N', 1, 'sp3', 2],
    ['N', 1, 'sp2', 1],
    ['N', 1, 'sp', 0],
    ['O', 2, 'sp3', 0],
    ['O', 1, 'sp3', 1],
    ['O', 1, 'sp2', 0],
    ['S', 4, 'sp3', 0],
    ['S', 3, 'sp3', 1],
    ['S', 3, 'sp2', 0],
    ['S', 2, 'sp3', 2],
    ['S', 2, 'sp2', 1],
    ['S', 1, 'sp2', 1],
    ['P', 4, 'sp3', 0],
    ['P', 3, 'sp3', 1],
    ['P', 3, 'sp2', 0],
    ['P', 2, 'sp3', 2],
    ['P', 2, 'sp2', 1],
    ['P', 2, 'sp', 0],
    ['P', 1, 'sp3', 3],
    ['B', 3, 'sp3', 0],
    ['B', 3, 'sp2', 0],
    ['B', 2, 'sp3', 1],
    ['B', 2, 'sp2', 1],
]

## [Atom 1, Atom 2, lower limit(bond length), upper limit(bond length), hybridization of atom 1]

bond_table = [
    [ "C", "C", 1.514, 1.588, 'sp3'], ## C-C
    [ "C", "C", 1.501, 1.513, 'sp2'], ## C=C
    [ "C", "C", 1.465, 1.473, 'sp'],  ## C-tb-C, tb = triple bond
    [ "C", "N", 1.487, 1.510, 'sp3'], ## C-N
    [ "C", "N", 1.278, 1.432, 'sp2'], ## C=N
    [ "C", "N", 1.134, 1.146, 'sp' ], ## C-tb-N
    [ "C", "O", 1.389, 1.50, 'sp3'],  ## C-O
    [ "C", "O", 1.180, 1.387, 'sp2'], ## C=O
    [ "C", "S", 1.744, 1.865, 'sp3'], ## C-S
    [ "C", "S", 1.710, 1.743, 'sp2'], ## C=S
    [ "C", "S", 1.629, 1.679, 'sp'],  ## C-tb-S
    [ "C", "P", 0.00, 2.00, 'sp3'],
    [ "C", "BR", 0.00, 3.00,'sp3'],
    [ "C", "CL", 0.00, 3.00, 'sp3'],
	[ "C", "Cl", 0.00, 3.00, 'sp3'],
    [ "C", "F", 0.00 , 3.00, 'sp3'],
    [ "C", "I", 0.00, 3.00, 'sp3'],
    [ "O", "C", 1.30, 3.00, 'sp3'],
    [ "O", "C", 0.00, 1.29, 'sp2'],
    [ "O", "S", 1.476, 3.00, 'sp3'],
    [ "O", "S", 0.00, 1.476, 'sp2'],
    [ "O", "P", 1.508, 3.00, 'sp3'],
    [ "O", "P", 0.00, 1.508, 'sp2'],
    [ "O", "N", 1.30, 3.00, 'sp3'],
    [ "O", "N", 0.00, 1.29, 'sp2'],
    [ "O", "O", 0.00, 3.00, 'sp3'],
    [ "N", "C", 1.468, 1.511, 'sp3'],
    [ "N", "C", 1.27, 1.432, 'sp2'],
    [ "N", "C", 1.13, 1.16, 'sp'],
    [ "N", "N", 1.300, 1.455, 'sp3'],
    [ "N", "N", 1.122, 1.299, 'sp2'],
    [ "N", "O", 1.24, 1.47, 'sp3'],
    [ "N", "O", 1.15, 1.239, 'sp2'],
    [ "N", "P", 1.625, 1.70, 'sp3'],
    [ "N", "P", 1.560, 1.61, 'sp2'],
    [ "N", "S", 1.59, 1.75, 'sp3'],
    [ "N", "S", 1.52, 1.58, 'sp2'],
    [ "N", "SI", 0.00, 3.00, 'sp3'],
	[ "N", "Si", 0.00, 3.00, 'sp3'],
    [ "P", "P", 2.00, 3.00, '???'],
    [ "P", "P", 0.00, 1.99, 'sp3'],
    [ "P", "N", 1.560, 1.61, 'sp3'],
    [ "P", "C", 0.00, 3.00, 'sp3'],
    [ "P", "O", 0.00, 3.00, 'sp3'],
    [ "P", "S", 0.00, 3.00, 'sp3'],
    [ "S", "N", 1.50, 1.6, 'sp'],
    [ "S", "N", 1.6, 3.00, 'sp2'],
    [ "S", "O", 1.50, 1.60, 'sp2'],
    [ "S", "O", 1.65, 3.00, 'sp'],
    [ "S", "N", 1.50, 3.00, 'sp2'],
    [ "S", "N", 0.00, 1.45, 'sp'],
    [ "S", "C", 0.00, 3.00, 'sp2']
]

## for use when you do not have to add hydrogens

hydrogen_table2 = [
    ['C', 4, 'sp3'],
    ['C', 3, 'sp2'],
    ['C', 2, 'sp'],
    ['N', 3, 'sp2'],   # Changed this from sp3 to sp2 by Salah
    ['N', 2, 'sp'],   # Changed this from sp2 to sp by Salah
    ['O', 2, 'sp3'],
    ['O', 1, 'sp2'],
    ['H', 1, 's'],
    ['S', 4, 'sp3'],
    ['P', 4, 'sp3']
]

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
	parser = argparse.ArgumentParser()
	parser.add_argument('whatever', action="store", help = 'file to parse')
	parser.add_argument('-p', action='store_false',dest='ideal',default =True, help = 'specifies file is Protein.pdb')
	parser.add_argument('-r', action='store_false',dest='reverse_order', default=True, help = 'print conformers starting from negative ions')
	options = parser.parse_args()
	
	pdb_file_name = fname.split('-')
	options.ligand = protonation_states.ligand_list[0][:3]
	options.charge = 1
	options.chain_identifier = 'A'
	
	
	i = 0
	j = 0
	pdb_file = pdb_file_name[0]+'.pdb'
	lines_model = open(pdb_file).readlines()
	list_of_files = []
	for line in lines_model:
		if lines_model[i].split()[0] == 'MODEL':
			start = True
			j += 1
		else: 
			start = False
		
		if start:
			file_model_name = pdb_file_name[0]+'_'+str(j)+'.pdb'
			list_of_files.append(file_model_name)
			fp =  open(file_model_name, 'w')
			fp.write(line)
		else:
			fp.write(line)
		i += 1
	
	k = 0
	for model_file in list_of_files:
		options.filename = model_file
		conf_chrage = protonation_states.ligand_list[k]
		conformers = mk_conformers(options,conf_chrage)
		pdb = Pdb(options)
		
		if options.ideal:
			add_hybrids(pdb)
		else:
			add_hydrogens(pdb)
		# con section
		write_con_section(options,tpl,pdb,conformers)
		k += 1
	
	
	
	# Atom Parameters
	write_atom_param_section(tpl,protonation_states,atoms_list,vdw_dict)
	
	# Atom Charges
	write_charges(tpl,protonation_states,atoms_list,numberofatoms)
	
	# BOND_ANG
	write_bond_angel(tpl,protonation_states,atoms_list,numberofatoms)
	
	# State Penalty
	write_extra(tpl,protonation_states,epik_State_Penalty)
	
	

	return


	
def main(fname_base,fname_pdb,fname, epik_State_Penalty_file): 

	lines = [line .strip() for line in open(fname).readlines()]
	if len(lines[1]) > 3:
		lines[1] = three_characters.upper()
		
	removeExt = fname_base.replace(".mol2", "")
	three_characters = removeExt[:3]
	paramDir = "param"
	if not os.path.exists(paramDir):
		os.makedirs(paramDir)
	#print 
	with open(os.path.join(paramDir, lines[1]+'.tpl'),'w') as tpl:
		write_tpl(tpl,fname,epik_State_Penalty_file)
	
	

if __name__ == "__main__":
	if len(sys.argv) < 2:
		print "Specify an input\nFor example python epik2tpl.py Bosutinib "
		sys.exit()
	fname_base = sys.argv[1]
	'''print len(fname_base.split('/'))
	if len(fname_base.split('/')) > 1:
		fname_merged = fname_base+'-merged.mol2'
		fname_state_penalties =  fname_base+'-state-penalties.out'
		fname_pdb = fname_base+'.pdb'
	else:'''	
	fname_merged = fname_base+'-merged.mol2'
	fname_state_penalties =  fname_base+'-state-penalties.out'
	fname_pdb = fname_base+'.pdb'
	
	#options = parser.parse_args()
	
	main(fname_base,fname_pdb,fname_merged,fname_state_penalties)
	
	
