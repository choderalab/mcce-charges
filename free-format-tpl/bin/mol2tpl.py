#!/usr/bin/python
import sys
from mapping import ORBITAL

def atom_in_residue(atom, residue):
    if residue.resName == atom.resName and \
            residue.resSeq == atom.resSeq:
        return True


def split_molecules(lines):
    molecules = []
    indexes = [i for i,x in enumerate(lines) if x == "@<TRIPOS>MOLECULE"]
    indexes.append(len(lines))
    if indexes:
        for i in range(1, len(indexes)):
            molecule = lines[indexes[i-1]:indexes[i]]
            molecules.append(molecule)

    return molecules


def same_xyz(v1, v2):
    if abs(v1[0] - v2[0]) < 0.001 and abs(v1[1] - v2[1]) < 0.001 and abs(v1[2] - v2[2]) < 0.001:
        return True
    else:
        return False


def identical_res(ref_residue, residue):
    identical_flag = True
    ref_atoms = ref_residue.atoms[:]
    atoms = residue.atoms[:]
    while atoms and ref_atoms:
        # pick the first one and try to get a match
        identical_found = False
        if atoms:
            for ref_atom in ref_atoms:
                if atoms[0].name == ref_atom.name and \
                        same_xyz(atoms[0].xyz, ref_atom.xyz) and \
                        abs(atoms[0].charge - ref_atom.charge) < 0.001:
                    # same atom and same xyz and charge
                    identical_found = True
                    atoms.pop(0)
                    ref_atoms.remove(ref_atom)
                    break

        if not identical_found:  # no reference atom is identical, a new conformer
            identical_flag = False
            break

    if (not atoms and not ref_atoms) and identical_flag: # one to one match all found and identical
        return True
    else:
        return False


def remove_identical_residue(residues):
    if residues:
        new_residues = [residues[0]]
        for residue in residues[1:]:
            is_new_residue = True
            for ref_residue in new_residues:
                if identical_res(ref_residue, residue):
                    is_new_residue = False
                    break

            if is_new_residue:
                new_residues.append(residue)

        return new_residues


class ATOM:
    def __init__(self, line):
        fields = line.split()
        self.serial = int(fields[0])
        self.name = fields[1]
        self.type = fields[5]
        self.resName = fields[7]
        self.resSeq = int(fields[6])
        self.xyz = (float(fields[2]), float(fields[3]), float(fields[4]))
        self.charge = float(fields[8])

        try:
            self.hybrid = ORBITAL[self.type]
        except:
            print "Error in processing this line:"
            print line
        self.connected = []

    def writeout(self):
        print "%4s %3s %4d %8.3f%8.3f%8.3f %8.3f" % \
              (self.name, self.resName, self.resSeq, self.xyz[0], self.xyz[1], self.xyz[2], self.charge)


class RESIDUE:
    def __init__(self, resname=None, resseq=None):
        self.resName = resname
        self.resSeq = resseq
        self.atoms = []


class PROTEIN:
    def __init__(self):
        self.residues = []

    def loadmol2(self, molecule):
        # read into atom lines into atom records
        atoms = []
        atom_records = False
        for line in molecule:
            if line.find("@<TRIPOS>ATOM") >= 0:
                atom_records = True
                continue
            elif line.find("@") >= 0:
                atom_records = False

            if atom_records:
                atom = ATOM(line)
                atoms.append(atom)

        # group into residues
        for atom in atoms:
            found = False
            for i_res in range(len(self.residues) - 1, -1, -1):
                if atom_in_residue(atom, self.residues[i_res]):
                    self.residues[i_res].atoms.append(atom)
                    found = True
                    break

            if not found:
                res = RESIDUE(resname=atom.resName, resseq=atom.resSeq)
                res.atoms.append(atom)
                self.residues.append(res)

        # Confirm serial number match index number
        for i in range(len(atoms)):
            if i != atoms[i].serial-1:
                atoms[i].writeout()
                print "Mismatched index: %d <-> %d" % (i, atoms[i].serial)

        # Construct bonds
        bond_records = False
        for line in molecule:
            if line.find("@<TRIPOS>BOND") >= 0:
                bond_records = True
                continue
            elif line.find("@") >= 0:
                bond_records = False

            if bond_records:
                fileds = line.split()
                i1 = int(fileds[1]) - 1
                i2 = int(fileds[2]) - 1
                atoms[i1].connected.append(atoms[i2])
                atoms[i2].connected.append(atoms[i1])


    def writetpl(self):
        for res in self.residues:
            print "Residue: %s %04d" % (res.resName, res.resSeq)
            for atom in res.atoms:
                connected = ", ".join(["\"%s\"" % x.name for x in atom.connected])
                line = "CONNECT, %s, \"%-4s\": %8.3f, radius, %5s, %s" % (atom.resName, atom.name, atom.charge,
                                                                        atom.hybrid,
                                                                         connected)
                print line




if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "Specify a PDB file and optional mol2 files with atom charges"
        sys.exit()

    fname = sys.argv[1]
    lines = [line .strip() for line in open(fname).readlines()]
    molecules = split_molecules(lines)

    proteins = []
    for molecule in molecules:
        prot = PROTEIN()
        prot.loadmol2(molecule)
        proteins.append(prot)
        #prot.writetpl()

    # Do tpl residue by residues
    tpl_residues = {}
    for prot in proteins:
        for residue in prot.residues:
            residue_id = "%s_%04d" % (residue.resName, residue.resSeq)
            if residue_id in tpl_residues:
                tpl_residues[residue_id].append(residue)
            else:
                tpl_residues[residue_id] = [residue]

    # process atoms based on residue
    for x in tpl_residues.keys():
        residues = tpl_residues[x]
        print len(residues),
        residues = remove_identical_residue(residues)
        print "-> %d" % len(residues)



