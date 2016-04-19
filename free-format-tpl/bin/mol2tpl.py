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

    for molecule in molecules:
        prot = PROTEIN()
        prot.loadmol2(molecule)
        prot.writetpl()
        print "===================="

