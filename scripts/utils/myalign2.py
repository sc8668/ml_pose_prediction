#!/usr/bin/env python
# -*- coding: utf-8 -*-

# =============================================================================
# align the structures using the structalign utility in Schrodinger 
# =============================================================================
import os, glob
from multiprocessing.dummy import Pool

def myoperation(t):
	com_files = glob.glob('%s/%s/*/*_prot/*_complex_prep.mae'%(os.getcwd(),t))
	cmd = 'cd %s &&'%t
	cmd += 'module load schrodinger &&'
	cmd += 'structalign %s'%' '.join(com_files)
	os.system(cmd)
	
	os.system("rm -rf %s/*/*_prot"%t)
	
	pdbs = [x for x in os.listdir('./%s'%t) if os.path.isdir('./%s/%s'%(t,x))]
	for pdb in pdbs:
		cmd = 'cd %s/%s &&'%(t, pdb)
		cmd += 'mkdir -p %s_prot &&'%pdb
		cmd += 'cd %s_prot &&'%pdb
		cmd += 'mv ../../rot-%s_complex_prep.mae %s_complex_prep.mae &&'%(pdb, pdb)
		cmd += 'module load schrodinger &&'
		cmd += 'run split_structure.py -many_files -k -m ligand %s_complex_prep.mae %s_complex_prep.mae &&'%(pdb, pdb)
		cmd += 'structconvert %s_complex_prep_receptor1.mae %s_p.pdb &&'%(pdb, pdb)
		cmd += 'structconvert %s_complex_prep_receptor1.mae %s_p.mol2 &&'%(pdb, pdb)
		cmd += 'structconvert %s_complex_prep_ligand1.mae %s_l.sdf &&'%(pdb, pdb)
		cmd += 'structconvert %s_complex_prep_ligand1.mae %s_l.mol2'%(pdb, pdb)
		os.system(cmd)

def main():
	targets = [x for x in os.listdir('.') if os.path.isdir(x)]
	pool = Pool(32)
	pool.map(myoperation, targets)
	pool.close()
	pool.join()


if __name__ == '__main__':
	main()	









def myoperation(t, pdb):
	cmd = 'cd %s/%s/%s_prot &&'%(t, pdb, pdb)
	cmd += 'module load schrodinger &&'
	cmd += 'run split_structure.py -many_files -k -m ligand %s_complex_prep.mae %s_complex_prep.mae &&'%(pdb, pdb)
	cmd += 'structconvert %s_complex_prep_receptor1.mae %s_p.pdb &&'%(pdb, pdb)
	cmd += 'structconvert %s_complex_prep_receptor1.mae %s_p.mol2 &&'%(pdb, pdb)
	cmd += 'structconvert %s_complex_prep_ligand1.mae %s_l.sdf &&'%(pdb, pdb)
	cmd += 'structconvert %s_complex_prep_ligand1.mae %s_l.mol2'%(pdb, pdb)
	os.system(cmd)

'''
names = [x for x in os.listdir('./293') if os.path.isdir('./293/%s'%(x))]
for name in names:
	os.system("cp 332/%s/%s_prot/%s_complex_prep.mae 332/"%(name,name,name))

total = ["293_%s_%s"%(x,y) for x in names for y in names if x!=y]
my=[]
x1=['1aid', '1b6j', '1b6k', '1b6l', '1d4k', '1d4l', '1mtr', '1z1h', '1kzk', '3aid', '3el1','3ekv', '3ekx','3mxd', '3mxe']
for x in x1:
	for y in x1:
		if x!=y:
			my.append("293_%s_%s"%(x,y))
		
final = [x for x in total if x not in my]

dfx = df.loc[names]
ww = [x for x in dfx.groupby('ligname') if len(x[-1])>1]
'''




