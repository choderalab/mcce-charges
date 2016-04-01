#!/usr/bin/python
import sys


def atom_in_residue(atom, residue):
    if residue.resName == atom.resName and \
            residue.chainID == atom.chainID and \
            residue.resSeq == atom.resSeq and \
            residue.iCode == atom.iCode:
        return True


class ATOM:
    def __init__(self, line):
        self.serial = int(line[6:11])
        self.name = line[12:16]
        self.altLoc = line[16]
        self.resName = line[17:20]
        self.chainID = line[21]
        self.resSeq = int(line[22:26])
        self.iCode = line[26]
        self.xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        self.connected = []

    def writeout(self):
        print "%4s %3s%c %c %4d%c %8.3f%8.3f%8.3f" % \
              (self.name, self.resName, self.altLoc, self.chainID, self.resSeq, self.iCode,
               self.xyz[0], self.xyz[1], self.xyz[2])


class RESIDUE:
    def __init__(self, resname=None, chainid=None, resseq=None, icode=None):
        self.resName = resname
        self.chainID = chainid
        self.resSeq = resseq
        self.iCode = icode
        self.atoms = []


class PROTEIN:
    def __init__(self):
        self.residues = []

    def loadpdb(self, fname):
        lines = open(fname).readlines()

        # read into atom lines into atom records
        atoms = []
        for line in lines:
            if line[:6] == "ATOM  " or line[:6] == "HETATM":
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
                res = RESIDUE(resname=atom.resName, chainid=atom.chainID, resseq=atom.resSeq, icode=atom.iCode)
                res.atoms.append(atom)
                self.residues.append(res)

    def writeout(self):
        for res in self.residues:
            print "Residue %s %c %4d%c" % (res.resName, res.chainID, res.resSeq, res.iCode)
            for atom in res.atoms:
                print "%s %8.3f%8.3f%8.3f" % (atom.name, atom.xyz[0], atom.xyz[1], atom.xyz[2])

    def connect(self):
        for res in self.residues:
            for atom in res.atoms:
                for res2 in self.residues:
                    for atom2 in res2.atoms:
                        if atom == atom2: pass
                        if bonded(atom, atom2):
                            atom.connected.append(atom2)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "Specify a PDB file"
        sys.exit()

    prot = PROTEIN()
    prot.loadpdb(sys.argv[1])
    # prot.writeout()
    prot.connect()

