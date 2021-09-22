#!/usr/bin/env python
# -*- coding: utf-8 -*-

# =============================================================================
# rescore with X-Score
# =============================================================================
import sys, os
import shutil
import subprocess
import multiprocessing
from multiprocessing import Manager
import pandas as pd


def write_file(output_file, outline):
	buffer = open(output_file, 'w')
	buffer.write(outline)
	buffer.close()


def xscore_calc(protname, ligname):
	'''for the docking poses'''
	cmdline = 'cd %s/%s_surflex/%s &&'%(protname, protname, ligname)
	cmdline += 'mkdir -p %s_xscore'%ligname
	os.system(cmdline)
	outline = '''FUNCTION    SCORE
RECEPTOR_PDB_FILE    ../../../%s_prot/%s_p.pdb
###prepare the input PDB file
FixPDB
REFERENCE_MOL2_FILE  ../../../%s_prot/%s_l.mol2 
#COFACTOR_MOL2_FILE  none 
LIGAND_MOL2_FILE     ../%s.mol2
###prepare the Mol2 files
FixMol2
#
OUTPUT_TABLE_FILE    ./xscore.table
OUTPUT_LOG_FILE      ./xscore.log
###
### how many top hits to extract from the LIGAND_MOL2_FILE?
###
NUMBER_OF_HITS       1 
HITS_DIRECTORY       ./xscore.mdb 
###
### want to include atomic binding scores in the resulting Mol2 files?
###
SHOW_ATOM_BIND_SCORE	NO		[YES/NO]
###
### set up scoring functions -----------------------------------------
###
APPLY_HPSCORE         YES             	[YES/NO]
	HPSCORE_CVDW  0.004 
	HPSCORE_CHB   0.054
	HPSCORE_CHP   0.009
	HPSCORE_CRT  -0.061
	HPSCORE_C0    3.441
APPLY_HMSCORE         YES             	[YES/NO]
	HMSCORE_CVDW  0.004
	HMSCORE_CHB   0.101
	HMSCORE_CHM   0.387
	HMSCORE_CRT  -0.097
	HMSCORE_C0    3.567
APPLY_HSSCORE         YES	  	[YES/NO]
	HSSCORE_CVDW  0.004
	HSSCORE_CHB   0.073
	HSSCORE_CHS   0.004
	HSSCORE_CRT  -0.090
	HSSCORE_C0    3.328
'''%(protname,protname,protname,protname,ligname)
	
	write_file('%s/%s_surflex/%s/%s_xscore/xscore.input'%(protname, protname, ligname, ligname), outline)
	##set the environment variable
	xscore_path='/home/shenchao/AI_based_SFs/software/XScore/xscore_v1.3/bin/xscore'
	cmd = 'export XSCORE_PARAMETER=/home/shenchao/AI_based_SFs/software/XScore/xscore_v1.3/parameter &&'
	cmd += '%s xscore.input'%xscore_path
	p = subprocess.Popen([cmd], shell=True, cwd="%s/%s_surflex/%s/%s_xscore"%(protname, protname, ligname, ligname))
	p.wait()
	
	try:
		lines = open('%s/%s_surflex/%s/%s_xscore/xscore.log'%(protname, protname, ligname, ligname), 'r').readlines()
		for line in lines:
			if line.startswith('Total'):
				VDW = line.split()[1].strip()
				HB = line.split()[2].strip()
				HP = line.split()[3].strip()
				HM = line.split()[4].strip()
				HS = line.split()[5].strip()
				RT = line.split()[6].strip()
				average_score = line.split()[-1].strip()
			if line.startswith('HPSCORE = '):
				HPSCORE = line.split('=')[-1].strip()
				HMSCORE = lines[lines.index(line)+1].split('=')[-1].strip()
				HSSCORE = lines[lines.index(line)+2].split('=')[-1].strip()
		results = [protname, ligname, average_score, HPSCORE, HMSCORE, HSSCORE, VDW, HB, HP, HM, HS, RT]
	except:
		results =  [protname, ligname, None, None, None, None, None, None, None, None, None, None]
	
	shutil.rmtree("%s/%s_surflex/%s/%s_xscore"%(protname, protname, ligname, ligname))
	return results


def xscore_calc0(protname):
	'''for the native poses'''
	cmdline = 'cd %s &&'%protname
	cmdline += 'mkdir -p %s_xscore'%protname
	os.system(cmdline)
	outline = '''FUNCTION    SCORE
RECEPTOR_PDB_FILE    ../%s_prot/%s_p.pdb
###prepare the input PDB file
FixPDB
REFERENCE_MOL2_FILE  ../%s_prot/%s_l.mol2 
#COFACTOR_MOL2_FILE  none 
LIGAND_MOL2_FILE     ../%s_prot/%s_l.mol2 
###prepare the Mol2 files
FixMol2
#
OUTPUT_TABLE_FILE    ./xscore.table
OUTPUT_LOG_FILE      ./xscore.log
###
### how many top hits to extract from the LIGAND_MOL2_FILE?
###
NUMBER_OF_HITS       1 
HITS_DIRECTORY       ./xscore.mdb 
###
### want to include atomic binding scores in the resulting Mol2 files?
###
SHOW_ATOM_BIND_SCORE	NO		[YES/NO]
###
### set up scoring functions -----------------------------------------
###
APPLY_HPSCORE         YES             	[YES/NO]
	HPSCORE_CVDW  0.004 
	HPSCORE_CHB   0.054
	HPSCORE_CHP   0.009
	HPSCORE_CRT  -0.061
	HPSCORE_C0    3.441
APPLY_HMSCORE         YES             	[YES/NO]
	HMSCORE_CVDW  0.004
	HMSCORE_CHB   0.101
	HMSCORE_CHM   0.387
	HMSCORE_CRT  -0.097
	HMSCORE_C0    3.567
APPLY_HSSCORE         YES	  	[YES/NO]
	HSSCORE_CVDW  0.004
	HSSCORE_CHB   0.073
	HSSCORE_CHS   0.004
	HSSCORE_CRT  -0.090
	HSSCORE_C0    3.328
'''%(protname,protname,protname,protname,protname,protname)
	
	write_file('%s//%s_xscore/xscore.input'%(protname, protname), outline)
	##set the environment variable
	xscore_path='/home/shenchao/AI_based_SFs/software/XScore/xscore_v1.3/bin/xscore'
	cmd = 'export XSCORE_PARAMETER=/home/shenchao/AI_based_SFs/software/XScore/xscore_v1.3/parameter &&'
	cmd += '%s xscore.input'%xscore_path
	p = subprocess.Popen([cmd], shell=True, cwd="%s/%s_xscore"%(protname, protname))
	p.wait()
	
	try:
		lines = open('%s/%s_xscore/xscore.log'%(protname, protname), 'r').readlines()
		for line in lines:
			if line.startswith('Total'):
				VDW = line.split()[1].strip()
				HB = line.split()[2].strip()
				HP = line.split()[3].strip()
				HM = line.split()[4].strip()
				HS = line.split()[5].strip()
				RT = line.split()[6].strip()
				average_score = line.split()[-1].strip()
			if line.startswith('HPSCORE = '):
				HPSCORE = line.split('=')[-1].strip()
				HMSCORE = lines[lines.index(line)+1].split('=')[-1].strip()
				HSSCORE = lines[lines.index(line)+2].split('=')[-1].strip()
		results = [protname, '1_000x', average_score, HPSCORE, HMSCORE, HSSCORE, VDW, HB, HP, HM, HS, RT]
	except:
		results =  [protname, '1_000x', None, None, None, None, None, None, None, None, None, None]
	
	#shutil.rmtree("%s/%s_xscore"%(protname, protname))
	return results



def xscore_calc2(name, i, return_dict):
	'''integrare the results'''
	lignames = [x for x in os.listdir('./%s/%s_surflex'%(name, name)) if (os.path.isdir('./%s/%s_surflex/%s'%(name, name, x)) and x!='1_000x')]
	results = []
	for ligname in lignames:
		result = xscore_calc(name, ligname)
		results.append(result)
	
	result0 = xscore_calc0(name)
	results.append(result0)	
	return_dict[i] = results



def main():
	names = [x for x in os.listdir('.') if os.path.isdir(x)]
	manger = Manager()
	return_dict = manger.dict()
	jobs = []
	pool = multiprocessing.Pool(32)
	for i, name in enumerate(names):
		p = pool.apply_async(xscore_calc2, args=(name, i, return_dict))
		jobs.append(p)
	
	pool.close()
	pool.join()

	df = pd.DataFrame(sum(return_dict.values(), []))
	df.columns = ['pdb_id', 'lig_id', 'average_score', 'HPSCORE', 'HMSCORE', 'HSSCORE', 'VDW', 'HB', 'HP', 'HM', 'HS', 'RT']
	df.sort_values(by=['pdb_id','lig_id'], inplace=True)
	df.reset_index(inplace=True)
	del df['index']
	df.to_csv("xscore_out.csv")




if __name__ == '__main__':
    main()














