#!/usr/bin/python
import sys
import os
import time
import os.path

from_ = '/home/salah/mcce-charges/mcce'
dirs = [d for d in os.listdir(from_) if os.path.isdir(os.path.join(from_, d))]

#dirs.remove('bin');
#dirs.remove('param04');

print dirs

for topDir in dirs:
	print topDir
	dirs_2 = [d for d in os.listdir(from_+'/'+str(topDir)) if os.path.isdir(os.path.join(from_+'/'+str(topDir), d))]
	sys_call = 'mkdir ' + topDir
	os.system(sys_call)
	for secondDir in dirs_2:
		dir_3 = [d for d in os.listdir(from_+'/'+str(topDir)+'/'+str(secondDir)) if os.path.isdir(os.path.join(from_+'/'+str(topDir)+ '/' +str(secondDir), d))]
		sys_call_2 = 'mkdir ' + topDir+'/'+secondDir
		os.system(sys_call_2)
		print sys_call_2
		for thirdDir in dir_3:
			dir_4 = [d for d in os.listdir(from_+'/'+str(topDir)+'/'+str(secondDir)+'/'+str(thirdDir)) if os.path.isdir(os.path.join(from_+'/'+str(topDir)+ '/' +str(secondDir)+'/'+str(thirdDir), d))]
			sys_call_4 = 'mkdir ' + topDir+'/'+secondDir+'/'+thirdDir
			os.system(sys_call_4)
			print thirdDir
			for fouthDir in dir_4:
				sys_call_5 = 'mkdir ' + topDir+'/'+secondDir+'/'+thirdDir+'/'+fouthDir
				os.system(sys_call_5)
				head3 = 'cp ' + from_ +'/'+ topDir+'/'+secondDir+'/'+thirdDir+'/'+fouthDir+'/head3.lst '+topDir+'/'+secondDir+'/'+thirdDir+'/'+fouthDir
				sum_charge = 'cp ' + from_ +'/'+ topDir+'/'+secondDir+'/'+thirdDir+'/'+fouthDir+'/sum_crg.out '+topDir+'/'+secondDir+'/'+thirdDir+'/'+fouthDir
				run_log = 'cp ' + from_ +'/'+ topDir+'/'+secondDir+'/'+thirdDir+'/'+fouthDir+'/run.log '+topDir+'/'+secondDir+'/'+thirdDir+'/'+fouthDir
				run_prm = 'cp ' + from_ +'/'+ topDir+'/'+secondDir+'/'+thirdDir+'/'+fouthDir+'/run.prm '+topDir+'/'+secondDir+'/'+thirdDir+'/'+fouthDir
				fort38 = 'cp ' + from_ +'/'+ topDir+'/'+secondDir+'/'+thirdDir+'/'+fouthDir+'/fort.38 '+topDir+'/'+secondDir+'/'+thirdDir+'/'+fouthDir
				
				os.system(head3)
				os.system(sum_charge)
				os.system(run_log)
				os.system(run_prm)
				os.system(fort38)
				
				
	
