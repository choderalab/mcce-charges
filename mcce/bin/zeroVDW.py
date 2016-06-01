#!/usr/bin/python
import sys
import re
import math
import os
import re
import fileinput




def replace_word(infile):
    old_word = '999.000'
    new_word = '000.000'
    if not os.path.isfile(infile):
        print ("Error on replace_word, not a regular file: "+infile)
        sys.exit(1)

    f1=open(infile,'r').read()
    f2=open(infile,'w')
    m=f1.replace(old_word,new_word)
    f2.write(m)

if __name__ == "__main__":
	if len(sys.argv) < 2:
		print "no "
		sys.exit()
	fname_base = sys.argv[1]
	replace_word(fname_base)
