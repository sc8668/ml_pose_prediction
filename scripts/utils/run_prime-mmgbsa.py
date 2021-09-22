#!/usr/bin/env python
# -*- coding: utf-8 -*-

# =============================================================================
# run Prime-MM/GBSA
# =============================================================================

import os, sys, glob, shutil
import multiprocessing
from multiprocessing.dummy import Pool
from multiprocessing import Manager
import pandas as pd
#from schrodinger import structure
	
	
def run_mmpbsa(protfile, ligfile, dist=0.0, ligflex= 'frozen', job_type='REAL_MIN', hostname='localhost', njobs=1):
	'''for the docking poses'''
	protname = os.path.basename(protfile).split('.')[0]
	ligname = os.path.basename(ligfile).split('.')[0]
	os.system('cd %s && mkdir -p %s_%s_prime-mmgbsa'%(os.path.dirname(ligfile), protname, ligname))
	
	#if not os.path.exists('%s/%s_%s_prime-mmgbsa/%s_%s_pv.maegz'%(os.path.dirname(ligfile), protname, ligname, protname, ligname)):
	#	get_pv_file(protfile, ligfile)
	
	cmdline = 'cd %s/%s_%s_prime-mmgbsa &&'%(os.path.dirname(ligfile), protname, ligname)
	##set the environment variable of Schrodinger.
	cmdline += 'module load schrodinger/2020-1 &&'
	cmdline += 'structcat -ipdb %s -isd %s -omae %s_%s_pv.maegz &&'%(protfile, ligfile, protname, ligname)
	if ligflex == 'frozen':
		cmdline += 'prime_mmgbsa -out_type PV -job_type REAL_MIN -report_prime_log yes -csv_output yes -OVERWRITE %s_%s_pv.maegz -HOST %s -NJOBS %s -prime_opt OPLS_VERSION=OPLS2005 -frozen -NOJOBID'%(protname, ligname, hostname, njobs)
	if ligflex == 'flexible':
		if float(dist) == 0.0:
			cmdline += 'prime_mmgbsa -out_type PV -job_type %s -report_prime_log yes -csv_output yes -OVERWRITE %s_%s_pv.maegz -HOST %s -NJOBS %s -prime_opt OPLS_VERSION=OPLS2005 -NOJOBID'%(job_type, protname, ligname, hostname, njobs)	
		if float(dist) > 0.0:
			cmdline += 'prime_mmgbsa -out_type COMPLEX -flexdist %.1f -job_type REAL_MIN -report_prime_log yes -csv_output yes -OVERWRITE %s_%s_pv.maegz -HOST %s -NJOBS %s -prime_opt OPLS_VERSION=OPLS2005 -NOJOBID'%(float(dist), protname, ligname, hostname, njobs)
	os.system(cmdline)
	
	

def run_mmpbsa2(name):
	'''for the native poses'''
	lignames = [x for x in os.listdir('./%s/%s_surflex'%(name, name)) if (os.path.isdir('./%s/%s_surflex/%s'%(name, name, x)) and x!='1_000x')]
	mypath = '/home/shenchao/AI_pose_filter/pdbbind2016'
	run_mmpbsa('%s/%s/%s_prot/%s_p.pdb'%(mypath,name,name,name), '%s/%s/%s_prot/%s_l.sdf'%(mypath,name,name,name))
	for ligname in lignames:
		run_mmpbsa('%s/%s/%s_prot/%s_p.pdb'%(mypath,name,name,name), '%s/%s/%s_surflex/%s/%s.sdf'%(mypath,name, name, ligname, ligname))
	

def integrate_results(name, i, return_dict):
	df_list = []
	df0 = pd.read_csv('%s/%s_prot/%s_p_%s_l_prime-mmgbsa/%s_p_%s_l-out.csv'%(name, name, name, name, name, name), header=0, index_col=None)
	df0.drop(labels=['title'], axis=1, inplace=True)
	df0['pdb_id'] = name
	df0['lig_id'] = '1_000x'
	df_list.append(df0)
	shutil.rmtree('%s/%s_prot/%s_p_%s_l_prime-mmgbsa'%(name, name, name, name))
	lignames = [x for x in os.listdir('./%s/%s_surflex'%(name, name)) if (os.path.isdir('./%s/%s_surflex/%s'%(name, name, x)) and x!='1_000x')]
	for ligname in lignames:
		df = pd.read_csv('%s/%s_surflex/%s/%s_p_%s_prime-mmgbsa/%s_p_%s-out.csv'%(name, name, ligname, name, ligname, name, ligname), header=0, index_col=None)
		df.drop(labels=['title'], axis=1, inplace=True)
		df['pdb_id'] = name
		df['lig_id'] = ligname
		df_list.append(df)
		shutil.rmtree('%s/%s_surflex/%s/%s_p_%s_prime-mmgbsa'%(name, name, ligname, name, ligname))
	return_dict[i] = df_list	


def integrate_results2(names):
	manger = Manager()
	return_dict = manger.dict()
	jobs = []
	pool = multiprocessing.Pool(32)
	for i, name in enumerate(names):
		p = pool.apply_async(integrate_results, args=(name, i, return_dict))
		jobs.append(p)
	
	pool.close()
	pool.join()	
	
	
	df0 = pd.concat(sum(return_dict.values(), []), axis=0)
	cols = list(df0)
	cols.insert(0, cols.pop(cols.index('lig_id')))
	cols.insert(0, cols.pop(cols.index('pdb_id')))
	df0 = df0[cols]
	df0.sort_values(by=['pdb_id','lig_id'], inplace=True)
	df0.reset_index(inplace=True)
	del df0['index']
	df0.to_csv("prime-mmgbsa_out.csv")
	os.system('bzip2 prime-mmgbsa_out.csv')
		
			

def main():
	names = [x for x in os.listdir('.') if os.path.isdir(x)]
	pool = Pool(32)
	pool.map(run_mmpbsa2, names)
	pool.close()
	pool.join()
	integrate_results2(names)



if __name__ == '__main__':
    main()



