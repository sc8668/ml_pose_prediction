#!/usr/bin/env python
# -*- coding: utf-8 -*-

# =============================================================================
# the demo to generate the ECIF fingerprint.
# =============================================================================

import sys, os
import subprocess
import multiprocessing
from multiprocessing import Manager
from multiprocessing.dummy import Pool
import pandas as pd


##the dir of ecif_calculate.py and PDB_Atom_Keys.csv
SCRIPT_dir = '~'

def ecif_calculate(pdbname, ligname):
	##set the envionment variables based on your own systems.
	cmd = 'module purge &&'
	cmd += 'export PYTHONPATH=/home/shenchao/python_module3.6/lib/python3.6/site-packages:/opt/anaconda3/5.2.0/lib/python3.6/site-packages &&'
	cmd += 'module load anaconda3/5.2.0 &&'
	if ligname == '1_000x':		
		cmd += 'python %s/ecif_calculate.py --pdb_atom_keys_file %s/PDB_Atom_Keys.csv -p ../../%s_prot/%s_p.pdb -l ../../%s_prot/%s_l.sdf --prot_name %s --lig_name %s -d ecif'%(SCRIPT_dir, SCRIPT_dir, pdbname, pdbname, pdbname, pdbname, pdbname, ligname)
		#cmd += 'python ./ecif_calculate.py --pdb_atom_keys_file %s/PDB_Atom_Keys.csv -p ../../%s_prot/%s_p.pdb -l ../../%s_prot/%s_l.sdf --prot_name %s --lig_name %s -d elem'%(SCRIPT_dir, SCRIPT_dir, pdbname, pdbname, pdbname, pdbname, pdbname, ligname)
	else:	
		cmd += 'python %s/ecif_calculate.py --pdb_atom_keys_file %s/PDB_Atom_Keys.csv -p ../../%s_prot/%s_p.pdb -l ./%s.sdf --prot_name %s --lig_name %s -d ecif'%(SCRIPT_dir, SCRIPT_dir, pdbname, pdbname, ligname, pdbname, ligname)
		#cmd += 'python ./ecif_calculate.py --pdb_atom_keys_file %s/PDB_Atom_Keys.csv -p ../../%s_prot/%s_p.pdb -l ./%s.sdf --prot_name %s --lig_name %s -d elem'%(SCRIPT_dir, SCRIPT_dir, pdbname, pdbname, ligname, pdbname, ligname)	
	p = subprocess.Popen([cmd], shell=True, cwd='%s/%s_surflex/%s'%(pdbname, pdbname, ligname))
	p.wait()	



def ecif_calculate2(pdbname):
	lignames = [x for x in os.listdir('./%s/%s_surflex'%(pdbname, pdbname)) if os.path.isdir('./%s/%s_surflex/%s'%(pdbname, pdbname, x))]
	for ligname in lignames:
		ecif_calculate(pdbname, ligname)


def integrate_results(name, i, return_dict, desc_type='ecif'):
	lignames = [x for x in os.listdir('./%s/%s_surflex'%(name, name)) if os.path.isdir('./%s/%s_surflex/%s'%(name, name, x))]
	df_list = []
	for ligname in lignames:
		df = pd.read_csv('%s/%s_surflex/%s/ecif_%s.csv'%(name, name, ligname, desc_type), header=0, index_col=0)
		os.remove('%s/%s_surflex/%s/ecif_%s.csv'%(name, name, ligname, desc_type))
		df_list.append(df)
	return_dict[i] = df_list
	

def integrate_results2(names, desc_type='ecif'):
	manger = Manager()
	return_dict = manger.dict()
	jobs = []
	pool = multiprocessing.Pool(32)
	for i, name in enumerate(names):
		p = pool.apply_async(integrate_results, args=(name, i, return_dict, desc_type))
		jobs.append(p)
	
	pool.close()
	pool.join()	
	
	
	df0 = pd.concat(sum(return_dict.values(), []), axis=0)
	df0.columns = ['pdb_id', 'lig_id', 'desc']
	df0.sort_values(by=['pdb_id','lig_id'], inplace=True)
	df0.reset_index(inplace=True)
	del df0['index']
	df0.to_csv("%s.csv"%desc_type)
	os.system('bzip2 %s.csv'%desc_type)



def main():
	names = [x for x in os.listdir('.') if os.path.isdir(x)]
	pool = Pool(32)
	pool.map(ecif_calculate2, names)
	pool.close()
	pool.join()
	
	integrate_results2(names, 'ecif')
	#integrate_results2(names, 'elem')



if __name__ == '__main__':
    main()


