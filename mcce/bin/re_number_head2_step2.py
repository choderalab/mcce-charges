#!/usr/bin/python
import sys
import re
import math
import commands
import os
import re

def deleteHH(head2_base, step2_base):
        print head2_base
	print step2_base
        head2_lines = open(head2_base).readlines()
	step2_lines = open(step2_base).readlines()
        
	for line in head2_lines:
                a = line
                b = head2_lines.next
                line = line.split()
		
                with open('new'+head2_base,'a') as tpl:
                        if line[1] == inhibitor:
                                tpl.write(a)
                                continue
                        if len(line[2]) != 4:
                                tpl.write(a)
                tpl.close()

if __name__ == "__main__":
        if len(sys.argv) < 3:
                print "need more input\nFormat:  python re_number_head2_step2.py head2.lst step2_out.pdb"
                sys.exit()
        head2_base = sys.argv[1]
        step2_base = sys.argv[2]
        deleteHH(head2_base,step2_base)

