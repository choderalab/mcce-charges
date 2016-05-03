#!/usr/bin/python
"""
This script takes the charge column from charged.mol2, put into mol2, print out the result.
"""
import sys

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "merge_charge.py HIS-HIS"
        print "   This will load HIS-HIS-epik-charged.mol2 and HIS-HIS-epik.mol2"
        sys.exit()

    fname_base = sys.argv[1]
    fname_charged = fname_base+"-epik-charged.mol2"
    fname_mol2 = fname_base+"-epik.mol2"

    lines_charged = open(fname_charged).readlines()
    lines_mol2 = open(fname_mol2).readlines()

    outputlines = []

    i_mol2 = 0
    i_charged_last = 0
    matched = False
    while i_mol2 < len(lines_mol2):
        # not matched, reset i_charged to last match
        i_charged = i_charged_last
        matched = False
        mol2_line = lines_mol2[i_mol2]
        if len(mol2_line) > 78:
            coord_mol2 = mol2_line[19:46]
            #print coord_mol2
            while i_charged < len(lines_charged):
                charged_line = lines_charged[i_charged]
                coord_charged = charged_line[19:46]
                if coord_mol2 == coord_charged:
                    newline = mol2_line[:70]+charged_line[70:77]+mol2_line[77:]
                    outputlines.append(newline)
                    i_charged_last = i_charged
                    i_charged += 1
                    i_mol2 += 1
                    matched = True
                    break
                else:
                    i_charged += 1

        if not matched:
            newline = mol2_line
            outputlines.append(newline)
            i_mol2 += 1

    for line in outputlines:
        print line,