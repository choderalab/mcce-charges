#!/usr/bin/python

## tpl_maker.py created by Denise M Kilburg for the Gunner Lab at CUNY-City College for the
## generation of ligand topology files for use in MCCE.

## tpl_maker.py requires a protein data bank file, ligand code to chain id, and charge.
## for example: tpl_maker.py 4Q6T.pdb "PO4 A" -3
## At the moment, the charge represents the most negative charge on an ion and it dictates
## the amount of conformers made. For example: in the above PO4 example, charge = -3 will
## create a PO401(neutral conformer), PO4-1, PO4-2,PO4-3 conformers. At the moment, the
## charge also specifies the number of hydrogens to be removed from the neutral conformer.
## At present, the program removes hydrogens from hydroxyl groups. If there are no
## hydroxyl groups it will remove hydrogens starting from H1, H2,...
## The add_hydrogens function adds hydrogens in 2 ways. If an atom is connected to 2 or more
## atoms, it looks at the angles that subtend them and assigns orbitals based on geometry: 
## sp3 == tetrahedral, sp2 == trigonal planar, sp == linear. It then adds the proper number of
## hydrogen atoms to fit the geometry. A Sulfur atom (S) for example
## will be assigned sp3 if it has tetrahedral geometry, not that it has sp3 hybridization.
## This needs to be corrected in the hydrogen_table, but for now this is how it is.
## If an atom only has one connection (terminal), it calculates bond distance and looks for 
## that bond distance in the bond_table and assigns a hybridization based on geometry.
## Send any questions or comments to:   dmkilburg@gmail.com

## To Chodera Lab: The formatting type functions are :
## header()
## natom()
## mk_iatom()
## mk_atomname()
## sect1_header()
## mk_proton()
## mk_pka()
## mk_electron()
## mk_EM()
## mk_RXN()
## mk_con_section()
## mk_atom_param_section()
## mk_charges()
## Please ignore the other functions. They are unfinished and ugly.



import sys
import math
import datetime


class atom(object):
    def __init__(self,name,id,e):
        self.name = name
        self.coords = []
        self.idnum = id
        self.connect = []
        self.element = e
        self.hybrid = ""
        
    def get_connects(self):
        return self.connect
    


class hydrogen(object):
    def __init__(self,name,connect,idnum):
        self.name = name
        self.connect = connect
        self.element = 'H'
        self.hybrid = 's'
        self.idnum = idnum

class conf(object):
    def __init__(self,name,charge):
        self.name = name
        self.charge = charge


######################################################

if len(sys.argv) < 4:
    print "Usage: tpl_maker.py pdbfile lig_code charge"
    sys.exit (1)


name = sys.argv[1]
ligand = sys.argv[2]
lig_pdb = ligand + ".pdb"
charge = int(sys.argv[3])
count = 0
lig_code = ligand[0:3]
order = 'r'
#order = ''
today = datetime.date.today()


################ Function Definitions ##################



def main_header(name,lig_code): ## Creates the commented out Main Header 
    output = open(lig_code+".tpl", 'w') 
    print >>output, "#"+ "-"*78 + "#"
    print >>output, "#" + " "*20 + "Topology File for:"
    print >>output, "#" + " "*28 +lig_code
    print >>output, "#" + " "*17 + "Extracted from: "+ name
    print >>output, "#" + " "*78 + "#"
    print >>output, "#" + " "*17 + "Created on: ", today
    print >>output, "#" + " "*78 + "#"
    print >>output, "#" + " "*17 + "Created with: tpl_maker.py"
    print >>output, "#" + "-"*78 + "#"
    output.close()
    return



def header(code,charge,order): ## Creates a Header (CONFLIST)
    output = open(code + ".tpl", 'a')
    conformers = []
    while charge < 0:
        con = code+str(charge)
        con = conf(con,charge)
        conformers.append(con)
        charge += 1
    con = code+"01"
    con = conf(con,charge)
    conformers.append(con)
    print >>output, "CONFLIST" + " " + code + " "*8 + code + "BK",
    if order == 'r':
        for i in reversed(conformers):
            print >>output, i.name,
    else:
        for i in conformers:
            print >>output, i.name,
    print >>output, code + "DM"    
    print >>output, '\n'
    output.close()
    return conformers

def mk_vectors(atom1,atom2,atom3):
    vect1 = []
    vect2 = []
    for i in range(3):
        vect1.append(float(atom2[i])-float(atom1[i]))
        vect2.append(float(atom3[i])-float(atom1[i]))
    return vect1,vect2


def dot_prod(vect1,vect2):
    dot_prod = 0
    for i in range(3):
        dot_prod = dot_prod + vect1[i]*vect2[i]
    return dot_prod

def magnitude(vector):
    mag = math.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)
    return mag


def find_angle(vect1,vect2):
        try:
            cos_ang = dot_prod(vect1,vect2)/(magnitude(vect1)*magnitude(vect2))
            angle = math.degrees(math.acos(cos_ang))
            return angle
        except:
            return None

## Bond distances from Allen, F.H.; Kennard, O. et al. (1987) Tables of 
## Bond Lengths determined by X-Ray and Neutron Diffraction. Part I . 
## Bond Lengths in Organic Compounds. J. CHEM. SOC. PERKIN TRANS

bond_table = [ ## [Atom 1, Atom 2, lower limit(bond length), upper limit(bond length), hybridization of atom 1]
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
def search_bond_table(a1,a2,distance):
    for row in bond_table:
        if row[0] == a1 and row[1] == a2 and row[2] <= distance and distance <= row[3]:
                return row[4]
    return 'unknown'

hydrogen_table = [ ## [Atom, number of atoms connected, hybridization, number of hydrogens to add]
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

def search_hydrogen_table(a1,connections, hybridization):
    for row in hydrogen_table:
        if row[0] == a1 and row[1] == connections and row[2] == hybridization:
            return row[3]
    return 0
   

def add_hydrogens(atom_list,con_dict): ## Work in progress
    j = 1
    hyd_list = []
    for atom in atom_list:
        con_list = atom.connect
        if len(con_list) >= 2:
            atom1 = atom.coords
            for elements in atom_list:
                if elements.idnum == con_list[0]:
                    atom2 = elements.coords
                elif elements.idnum == con_list[1]:
                    atom3 = elements.coords
                else:
                    pass
            (vector1, vector2) = mk_vectors(atom1,atom2,atom3)
            angle = find_angle(vector1,vector2)
            atom.hybrid = geometry(angle)
            num_hydrogens = search_hydrogen_table(atom.element,len(con_list),atom.hybrid)
            if num_hydrogens == 0:
                pass
            else:
                for i in range(num_hydrogens):
                    atom.connect.append(j)
                    (con_dict,hyd_list,j) = add_hydros(j,con_dict,atom.name,hyd_list)
        elif len(con_list) == 1:
            for elements in atom_list:
                if elements.idnum == con_list[0]:
                    bond = bond_length(atom.coords,elements.coords)
                    orbital = search_bond_table(atom.element,elements.element,bond)
                    atom.hybrid = orbital
                    num_hydrogens = search_hydrogen_table(atom.element,len(con_list),atom.hybrid)
                    if num_hydrogens == 0:
                        pass
                    else:
                        for i in range(num_hydrogens):
                            atom.connect.append(j)
                            (con_dict,hyd_list,j) = add_hydros(j,con_dict,atom.name,hyd_list)
                else:
                    pass
        else:
            pass
    return con_dict, hyd_list



def bond_length(atom1,atom2):
    vector = []
    for i in range(3):
        vector.append(float(atom1[i])-float(atom2[i]))
    length = magnitude(vector)
    return length

def add_hydros(j,con_dict,name,hyd_list): ## Work in progress
    
    con_dict[j] = 'H'+str(j)
    hyd = 'H'+str(j)
    hyd = hydrogen(hyd,name,j)
    hyd_list.append(hyd)
    j+= 1
    return con_dict,hyd_list,j

def geometry(angle): ## Work in progress
    if angle > 106 and angle < 115:
        return 'sp3'
    if angle < 183 and angle > 177:
        return 'sp'
    if angle < 125 and angle > 115:
        return 'sp2'
    else:
        return 'unknown'

def natom(conformers,atom_list,code,charge,hyd_list,order):
    max_atoms = len(atom_list) + len(hyd_list)
    min_atoms = max_atoms + charge
    template = '{0:9}{1:10}{2:1}\n'
    output = open(code + ".tpl", 'a')
    output.write(template.format('NATOM',code+"BK", '0'))
    if order == 'r':
        for conformer in reversed(conformers):
            output.write(template.format('NATOM',conformer.name, str(max_atoms)))
            max_atoms -= 1
    else:
        for conformer in conformers:
            output.write(template.format('NATOM',conformer.name, str(min_atoms)))
            min_atoms += 1
    output.write(template.format('NATOM',code+"DM", '0'))
    output.write('\n')
    output.close
    return



def mk_lig_file(name, lig_pdb): ## grep-like function. Makes a file (ligand.pdb)
    file = open(name)
    output = open(lig_pdb, 'w')
    for line in file:
        if "HETATM" in line:
            if ligand in line:
                print >> output, line,
    output.close()
    file.close()
    return


def mk_atoms(lig_pdb): 
    file = open(lig_pdb, 'r')
    atom_list = []
    for line in file:
        column = line.split()
        column0 = column[0]
        if len(column0) > 6:
            column[1] = atom(column[1],column0[6:],column[10])
            atom_list.append(column[1])
        else:
            column[2]  = atom(column[2],column[1],column[11])
            atom_list.append(column[2])
    file.close()
    return atom_list

def mk_connex(lig_code,atom_list):
    file = open(lig_code+"_con.txt", 'r')
    i = 0
    for line in file:
        column = line.split()
        for j in range(1,len(column)):
            atom_list[i].connect.append(column[j])
        i += 1
    file.close()
    return

def mk_coords(lig_pdb,atom_list): 
    file = open(lig_pdb, 'r')     
    i = 0
    for line in file:
        column = line.split()
        for j in range(6,9):
            atom_list[i].coords.append(column[j])
        i += 1
    file.close()
    return


def mk_iatom(lig_code,atom_list,charge,hyd_list,conformers,order):
    output = open(lig_code+".tpl", 'a')
    template = '{0:9}{1:7}{2:4}{3:4}\n'
    if order == 'r':
        conformers = reversed(conformers)
    else:
        pass
    for conformer in conformers:
        count = 0
        for atom in atom_list:
            output.write(template.format('IATOM',conformer.name, atom.name, str(count)))
            count += 1
            j = conformer.charge
        for hyd in hyd_list:
            if j < 0:
                j += 1
            else:
                output.write(template.format('IATOM',conformer.name, hyd.name, str(count)))
                count += 1
        output.write('\n')
    output.close()
    return


def mk_atomname(lig_code,atom_list,charge,hyd_list,conformers,order):
    output = open(lig_code+".tpl", 'a')
    template = '{0:9}{1:8}{2:>2}  {3:5}\n'
    if order == 'r':
        conformers = reversed(conformers)
    else:
        pass
    for conformer in conformers:
        count = 0
        for atom in atom_list:
            output.write(template.format("ATOMNAME",conformer.name,str(count),atom.name))
            count += 1
        j = conformer.charge
        for hyd in hyd_list:
            if j < 0:
                j += 1
            else:
                output.write(template.format("ATOMNAME",conformer.name,str(count),hyd.name))
                count += 1
        output.write('\n')
    output.close()    
    return

def get_con(name,lig_pdb,code): ## writes the _con.txt file that is needed 
    output = open(code + "_con.txt", 'w') ## for the get_connects function
    connect = open(lig_pdb)
    for line in connect:
        column = line.split()
        res_num = column[1]
        with open(name, 'r') as source:
            for row in source:
                if ("CONECT "+res_num) in row:
                    rcolumn = row.split()
             
                    r0column = rcolumn[0]
                    print >>output, r0column[6:] + " ".join(rcolumn[1:])
                elif ("CONECT"+res_num) in row:
                    rcolumn = row.split()
                    print >>output, " ".join(rcolumn[1:])
    connect.close()
    output.close()
    return


def create_con_dict(lig_pdb): ## Creates a dictionary that assigns an
    con_dict = {}             ## atom name to a res id number
    with open(lig_pdb, 'r') as f:
        for line in f:
            column = line.split()
            number = column[1]
            atom = column[2]
            con_dict[number] = atom

    return con_dict


def sect1_header(lig_code): 
    output =  open(lig_code + ".tpl", 'a') 
    output.write("# 1. Basic Conformer Information:\n")
    output.write("# Number of protons and electrons, pKa, Em, and Reaction Field Energy (RXN)\n")
    output.write('\n')
    output.close()
    return



def mk_proton(lig_code,conformers,order):
    output = open(lig_code+".tpl", 'a')
    output.write("# PROTON SECTION: PROTON means charge:\n")
    template = '{0:9}{1:11}{2:5}\n'
    if order == 'r':
        conformers = reversed(conformers)
    else:
        pass
    for conformer in conformers:
        output.write(template.format("PROTON",conformer.name,str(conformer.charge)))
    
    output.write('\n')

    output.close()
    return
    


def mk_pka(lig_code,conformers,order): 
    output = open(lig_code+".tpl", 'a')
    output.write("# Solution pKa Section: pKa data from CRC Handbook of Chemistry and Physics\n")
    template = '{0:9}{1:11}{2:5}\n'
    if order == 'r':
        conformers = reversed(conformers)
    else:
        pass
    for conformer in conformers:
        output.write(template.format("PKA", conformer.name, "0.0"))
    output.write(template.format("PKA", lig_code+"DM", "0.0"))
    output.write('\n')
    output.close()
    return


def mk_electron(lig_code,conformers,order):
    output = open(lig_code+".tpl", 'a')
    output.write("#ELECTRON SECTION:\n")
    template = '{0:9}{1:11}{2:5}\n'
    if order == 'r':
        conformers = reversed(conformers)
    else:
        pass
    for conformer in conformers:
        output.write(template.format("ELECTRON",conformer.name,"0.0"))
    
    output.write('\n')
    output.close()
    return

def mk_EM(lig_code,conformers,order): 
    output = open(lig_code+".tpl", 'a')
    template = '{0:9}{1:11}{2:5}\n'
    output.write("# EM SECTION:\n")
    if order == 'r':
        conformers = reversed(conformers)
    else:
        pass
    for conformer in conformers:
        output.write(template.format("EM",conformer.name,"0.0"))
    
    output.write('\n')
    
    output.close()
    return

def mk_RXN(lig_code): 
    output = open(lig_code+".tpl", 'a')
    template = '{0:9}{1:11}{2:5}\n'
    output.write("# REACTION FIELD ENERGY SECTION:\n")
    if order == 'r':
        conformers = reversed(conformers)
    else:
        pass
    for conformer in conformers:
        output.write(template.format("RXN",conformer.name,"0.0"))
    
    output.write('\n')
    output.close()
    return

def con_header(output): 

    string = "ires" + " " + "conn" + " "
    print >>output, "#ONNECT" + "   " + "conf" + " " + "atom" + "  " + "orbital" + "  " +string*4
    print >>output, "#ONNECT" +" " + "|-----|----|---------|" + "----|"*10

    return

def mk_con_section(atom_list,lig_code,con_dict,conformers,order):
    output = open(lig_code+".tpl",'a')
    if order == 'r':
        conformers = reversed(conformers)
    else:
        pass
    for conformer in conformers:
        print >>output, "#  "+ conformer.name
        con_header(output)
        template1 = '{0:9}{1:5} {2:^4} {3:^10}'
        template2 = '{0:^4} {1:^4} '
        template3 = '{0:9}{1:5} {2:^4} {3:^10}{4:^4} {5:^4}\n'
        j = conformer.charge
        for atom in atom_list:
            con_list = atom.get_connects()
            output.write(template1.format("CONNECT",conformer.name,atom.name,atom.hybrid))
            if j < 0:
                if atom.element == 'O' and atom.hybrid == 'sp3':
                    for connects in con_list[:-1]:
                        output.write(template2.format("0",con_dict[connects]))
                    j += 1
                else:
                    for connects in con_list:
                        output.write(template2.format("0",con_dict[connects]))
            else:
                for connects in con_list:
                    output.write(template2.format("0",con_dict[connects]))
            output.write('\n')
        j = conformer.charge
        for hyd in hyd_list:
            if j < 0:
                j += 1
            else:
                output.write(template3.format("CONNECT",conformer.name,hyd.name,hyd.hybrid,"0",hyd.connect))

        output.write('\n')
    output.close()
    return

def mk_atom_param_section(lig_code,charge,atom_list,vdw_dict,hyd_list,conformers,order):
    output = open(lig_code+".tpl", 'a')
    output.write("# Atom Parameters:\n")
    output.write("# Van Der Waals Radii. See source for reference\n")
    template = '{0:9}{1:7}{2:3} {3:7}\n'
    if order == 'r':
        conformers = reversed(conformers)
    else:
        pass
    for conformer in conformers:
        for atom in atom_list:
            output.write(template.format("RADIUS",conformer.name,atom.name,str(vdw_dict[atom.element])))
        j = conformer.charge    
        for hyd in hyd_list:
            if j < 0:
                j += 1
            else:
                output.write(template.format("RADIUS",conformer.name,hyd.name,str(vdw_dict['H'])))
        output.write('\n')
    output.write('\n')
    output.close()
    return

def mk_charges(lig_code,charge,atom_list,hyd_list,conformers,order): ## Work in progress
    output = open(lig_code+".tpl", 'a')
    template = '{0:9}{1:7}{2:3} {3:7}\n'
    if order == 'r':
        conformers = reversed(conformers)
    else:
        pass
    for conformer in conformers:
        i = 0
        try:
            charges = parse_g09(conformer.name+'chrg.log')
            for atom in atom_list:
                output.write(template.format("CHARGE",conformer.name,atom.name,charges[i][2]))
                i += 1
            j = conformer.charge
            for hyd in hyd_list:
                if j < 0:
                    j += 1
                else:
                    output.write(template.format("CHARGE",conformer.name,hyd.name,charges[i][2]))
                    i += 1
        except:
            for atom in atom_list:
                output.write(template.format("CHARGE",conformer.name,atom.name,'0.0'))
            j = conformer.charge
            for hyd in hyd_list:
                if j < 0:
                    j +=1
                else:
                    output.write(template.format("CHARGE",conformer.name,hyd.name,'0.0'))
        output.write('\n')
    output.write('\n')
    output.close()
    return

import re
 
def parse_g09(filename):
    found_charges=0
    charges=[]
    rex = re.compile("\s+([1-9][0-9]*)\s+([A-Za-z]{1,2})\s+(-{0,1}[0-9]+.[0-9]+)")
    for line in open(filename):
        if found_charges:
            if '------------' in line:
                return charges
            m = re.match(rex,line)
            if m != None:
                charges.append(m.groups())
            elif len(charges) > 0:
                print "Expected to match line: " + line
        elif 'Fitting point charges to electrostatic potential' in line:
            found_charges=1

    print 'Did not expect to reach the end of file ' + filename
    sys.exit(1)
   

#################  Function  ###########################################
######################## Calls ####################################

vdw_dict = {'H':'1.20','He':'1.40','C':'1.70','N':1.55,'O':1.52,'F':1.47,'NE':1.54,'SI':2.10,'P':'1.80','S':'1.80','CL':1.75,'AR':1.88,'AS':1.85,'SE':'1.90','BR':1.85,'KR':2.02,'TE':2.06,'I':1.98,'XE':2.16,'ZN':1.39,'CU':'1.40','HG':1.55,'CD':1.58,'NI':1.63,'PD':1.63,'AU':1.66,'AG':1.72,'MG':1.73,'PT':1.75,'LI':1.82,'U':1.86,'GA':1.87,'PB':2.02,'SN':2.17,'NA':2.27,'K':2.75}
##vdw_dict contains vdw radius of atoms as written in Bondi,A. van Der Waals Volumes and Radii,
##J.Phys.Chem, 63, 3, 1964 and http://periodic.lanl.gov/index.shtml (Los Alamos National Lab Periodic Table)


mk_lig_file(name,lig_pdb)

main_header(name,lig_code)

conformers = header(lig_code,charge,order)

atom_list = mk_atoms(lig_pdb)

get_con(name,lig_pdb,lig_code)

con_dict = create_con_dict(lig_pdb)

mk_connex(lig_code,atom_list)

mk_coords(lig_pdb,atom_list)

(con_dict, hyd_list) = add_hydrogens(atom_list,con_dict)
## The functions below this line can be commented out
## if you do not want a certain section to be printed. However, you cannot
## comment out the functions above, as all the functions below are dependent on
## them.

natom(conformers,atom_list,lig_code,charge,hyd_list,order)

mk_iatom(lig_code,atom_list,charge,hyd_list,conformers,order)

mk_atomname(lig_code,atom_list,charge,hyd_list,conformers,order)

sect1_header(lig_code)

mk_proton(lig_code,conformers,order)

mk_pka(lig_code,conformers,order)

mk_electron(lig_code,conformers,order)

mk_EM(lig_code,conformers,order)

mk_RXN(lig_code)

mk_con_section(atom_list,lig_code,con_dict,conformers,order)

mk_atom_param_section(lig_code,charge,atom_list,vdw_dict,hyd_list,conformers,order)

mk_charges(lig_code,charge,atom_list,hyd_list,conformers,order)



