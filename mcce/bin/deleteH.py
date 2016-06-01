#!/usr/bin/python
import sys
import re
import math
import datetime
import commands
import os
import re

def deleteHH(fname_base, inhibitor):
	print fname_base
	lines_ = open(fname_base).readlines()
	for line in lines_:
		a = line
		#print a
		line = line.split()
		with open('new'+fname_base,'a') as tpl:
			if line[3] == inhibitor:
				tpl.write(a)
				continue
			if line[2][:1] != 'H':
				tpl.write(a)		
		tpl.close()

if __name__ == "__main__":
	if len(sys.argv) < 3:
		print "need more input\nFormat:  python deleteH.py 4LMN_fixed_ph7.4.pdb EUI "
		sys.exit()
	fname_base = sys.argv[1]
	inhibitor = sys.argv[2]
	deleteHH(fname_base,inhibitor)
