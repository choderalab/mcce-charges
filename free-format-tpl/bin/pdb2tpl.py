#!/usr/bin/python
import sys
from param import COVALENT_RADIUS
from param import ATOM_MAPPING
from geom import *

UPPER109 = (109.4712+5.0)/180.0*math.pi
LOWER109 = (109.4712-5.0)/180.0*math.pi
UPPER120 = (120+5.0)/180.0*math.pi
LOWER120 = (120-5.0)/180.0*math.pi
UPPER180 = (180+5.0)/180.0*math.pi
LOWER180 = (180-5.0)/180.0*math.pi


def atom_in_residue(atom, residue):
    if residue.resName == atom.resName and \
            residue.chainID == atom.chainID and \
            residue.resSeq == atom.resSeq and \
            residue.iCode == atom.iCode:
        return True


def bonded(atom1, atom2):
    r1 = COVALENT_RADIUS[atom1.name]
    r2 = COVALENT_RADIUS[atom2.name]
    bd = (r1+r2) * 1.1   # 1.1 for room of error
    bd2 = bd*bd

    dx = atom1.xyz.x - atom2.xyz.x
    if abs(dx) > bd:
        return False
    dy = atom1.xyz.y - atom2.xyz.y
    if abs(dy) > bd:
        return False
    dz = atom1.xyz.z - atom2.xyz.z
    if abs(dz) > bd:
        return False
    dd = dx*dx + dy*dy + dz*dz
    if dd > bd2:
        return False

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
        self.xyz = VECTOR(float(line[30:38]), float(line[38:46]), float(line[46:54]))
        self.connected = []
        self.hybrid = "unk"
        if self.name in ATOM_MAPPING:
            self.name = ATOM_MAPPING[self.name]

    def writeout(self):
        print "%4s %3s%c %c %4d%c %8.3f%8.3f%8.3f" % \
              (self.name, self.resName, self.altLoc, self.chainID, self.resSeq, self.iCode,
               self.xyz.x, self.xyz.y, self.xyz.z)


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

        self.connect()
        self.hybrid()

    def hybrid(self):
        """
        Bonds = 6 -> sp3d2
        Bonds = 5 -> undecided
        Bonds = 4, -> bond angle = 109 -> sp3
                   |
                    -> undecided
        Bonds = 4, -> bond angle = 109 -> sp3
                   |
                    -> undecided
        Bonds = 3, -> bond angle = 109 -> sp3
                   |
                    -> bond angle = 120 -> sp2
        Bonds = 2, -> bond angle = 109 -> sp3
                   |
                    -> bond angle = 120 -> sp2
                   |
                    -> bond angle = 180 -> sp
        Bonds = 1, -> H -> s
                   |
                    -> undecided
        """
        for res in self.residues:
            for atom in res.atoms:
                bonds = len(atom.connected)
                if bonds == 0:
                    pass
                elif bonds == 1:
                    if atom.name[0] == "H" or atom.name[:2] == " H":
                        atom.hybrid = "s"
                    elif atom.name[:2] == " O":
                        atom.hybrid = "sp2"
                else:
                    v1 = vector_vminusv(atom.connected[0].xyz, atom.xyz)
                    v2 = vector_vminusv(atom.connected[1].xyz, atom.xyz)
                    alpha = avv(v1, v2)

                    if bonds == 6:
                        atom.hybrid = "sp3d2"
                    elif bonds == 4:
                        if UPPER109 > alpha > LOWER109: atom.hybrid = "sp3"
                        else: atom.hybrid = "%.1f" % alpha
                    elif bonds == 3:
                        if UPPER109 > alpha > LOWER109:
                            atom.hybrid = "sp3"
                        elif UPPER120 > alpha > LOWER120:
                            atom.hybrid = "sp2"
                        else:
                            atom.hybrid = "%.1f" % alpha
                    elif bonds == 2:
                        if UPPER109 > alpha > LOWER109:
                            atom.hybrid = "sp3"
                        elif UPPER120 > alpha > LOWER120:
                            atom.hybrid = "sp2"
                        elif UPPER180 > alpha > LOWER180:
                            atom.hybrid = "sp"
                        else:
                            atom.hybrid = "%.1f" % alpha

        return


    def writetpl(self):
        for res in self.residues:
            print "Residue: %s %c %04d" % (res.resName, res.chainID, res.resSeq)
            for atom in res.atoms:
                connected = ", ".join(["\"%s\"" % x.name for x in atom.connected])
                line = "CONNECT, %s, \"%s\": charge, radius, %5s, %s" % (atom.resName, atom.name, atom.hybrid,
                                                                         connected)
                print line


    def connect(self):
        for res in self.residues:
            for atom in res.atoms:
                for res2 in self.residues:
                    for atom2 in res2.atoms:
                        if atom != atom2 and bonded(atom, atom2):
                            atom.connected.append(atom2)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "Specify a PDB file"
        sys.exit()

    prot = PROTEIN()
    prot.loadpdb(sys.argv[1])
    prot.writetpl()

