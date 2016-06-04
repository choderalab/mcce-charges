#!/usr/bin/python
import sys
import os
import time
import os.path

dirs = [d for d in os.listdir('.') if os.path.isdir(os.path.join('.', d))]
for the_dir in dirs:
	print '***************************************************'
	if the_dir == 'param' or the_dir == 'bin':
		continue
	print 'Running merge_charge.py on ' + the_dir
	sys_call = 'python bin/merge_charge.py '+the_dir+'/'+the_dir+' > '+the_dir+'/'+the_dir+'-merged.mol2'
	os.system(sys_call)
	time.sleep( 0.2 )
	
	
	print 'Running epik2tpl_1.py on ' + the_dir
	sys_call = 'python bin/epik2tpl_1.py '+the_dir+'/'+the_dir
	os.system(sys_call)
	print '\n\n'

	
