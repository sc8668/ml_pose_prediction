#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ========================================================================
# a script to generate features
# ========================================================================

import os, sys, re, csv
import gzip, bz2
import argparse
import numpy as np
import pandas as pd
import subprocess
import multiprocessing
from multiprocessing import Manager
import uuid

MGLTOOLS_DIR="/opt/mgltools/1.5.6"
BABEL_DIR="/opt/openbabel/3.1.0/bin"
VINA_EXEC="/opt/vina/1.1.2/bin/vina"
SCRIPT_DIR="~/data/scripts"


uuid_str = uuid.uuid4().hex

def UserInput():
	p = argparse.ArgumentParser(description='Command Line Arguments')
	p.add_argument('-p','--prot', required=True,
				help='Input protein file (.pdb)')
	p.add_argument('-l','--lig', required=True,
				help='Input ligand docking poses (.mol2 or .sdf)  * gzip/bzip2 accepted,\
				Now only a single ligand with multiple docking poses is supported.\
				The poses should be ranked by their docking scores.')
	p.add_argument('-pn', '--protname', default=None, help="The name of the protein.")
	p.add_argument('-f', '--feats', default = 'all', choices=['vina','ecif','all'], help="The features to generate.")
	p.add_argument('-o','--outname', default = None, help='The prefix of the output files.')
	p.add_argument('-bo','--bzip_output', action="store_true", default = False, help='whether to compress the output file using "bzip2"')
	args = p.parse_args()
	
	args.prot = os.path.abspath(args.prot)
	args.lig = os.path.abspath(args.lig)
	
	if args.protname is None:
		args.protname = os.path.basename(args.prot).split('.')[0]
	if args.outname is None:
		args.outname = os.path.basename(args.prot).split('.')[0]		
	return args


def write_file(output_file, outline):
    buffer = open(output_file, 'w')
    buffer.write(outline)
    buffer.close()


def file_handle(file_name):
	'''
	Handle gzip and bzip2 file if the extension is right. otherwise, just open
	'''
	if re.search(r'.gz$', file_name):
		with gzip.open(file_name, 'rb') as f:
			contents = re.sub('\r\n', '\n', f.read().decode('utf-8'), count=0)			
	elif re.search(r'.bz2$', file_name):
		with bz2.open(file_name, 'rb') as f:
			contents = re.sub('\r\n', '\n', f.read().decode('utf-8'), count=0)
	else:
		contents = open(file_name, 'r').read()
	return contents


def prot_prep(prot, protname):
	'''prepare the protein'''
	cmd = "%s/bin/pythonsh "%MGLTOOLS_DIR
	cmd += "%s/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py "%MGLTOOLS_DIR
	cmd += "-r %s -o %s.pdbqt  -A checkhydrogens -U nphs_lps_waters"%(prot, protname)
	os.system(cmd)


def ligs_split(lig):
	'''split the ligands'''
	contents = file_handle(lig)
	os.system("mkdir -p temp_%s"%uuid_str)
	if not re.search(r'.sdf(.*)$|.mol2(.*)$', lig):
		print("Only the ligands with .mol2 or .sdf are supported !")
		sys.exit(1)
	elif re.search(r'.sdf(.*)$', lig):
		if re.search(r'.bz2$', lig):
			os.system("bunzip2 %s"%lig)
			cmd = "%s/obabel -isd %s -omol2 | sed '/@<TRIPOS>UNITY_ATOM_ATTR/,/@<TRIPOS>BOND/c@<TRIPOS>BOND'> %s.mol2"%(BABEL_DIR, lig.replace(".bz2",""), os.path.basename(lig).split('.')[0])
			
		else:
			cmd = "%s/obabel -isd %s -omol2 | sed '/@<TRIPOS>UNITY_ATOM_ATTR/,/@<TRIPOS>BOND/c@<TRIPOS>BOND'> %s.mol2"%(BABEL_DIR, lig, os.path.basename(lig).split('.')[0])
		os.system(cmd)
		contents = open("%s.mol2"%os.path.basename(lig).split('.')[0], 'r').read().split('@<TRIPOS>MOLECULE\n')
	
	else:
		if re.search(r'.gz$', lig):
			with gzip.open(lig, 'rb') as f:
				contents = re.sub('\r\n', '\n', f.read().decode('utf-8'), count=0).split('@<TRIPOS>MOLECULE\n')
			
		#elif re.search(r'.bz2$', lig):
			with bz2.open(lig, 'rb') as f:
				contents = re.sub('\r\n', '\n', f.read().decode('utf-8'), count=0).split('@<TRIPOS>MOLECULE\n')
		else:
			contents = open(lig, 'r').read().split('@<TRIPOS>MOLECULE\n')
	for i, c in enumerate(contents[1:]):
		os.system("mkdir -p temp_%s/lig_%s"%(uuid_str, i))
		write_file('temp_%s/lig_%s/lig_%s.mol2'%(uuid_str,i,i), '@<TRIPOS>MOLECULE\n'+c)


def lig_prep(lig, ligname, ligdir='.'):
	'''prepare the ligand'''
	
	cmd = "%s/bin/pythonsh "%MGLTOOLS_DIR
	cmd += "%s/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py "%MGLTOOLS_DIR
	cmd += "-l %s -o %s.pdbqt -A checkhydrogens"%(lig, ligname)
	p = subprocess.Popen([cmd], shell=True, cwd=ligdir)
	p.wait()
	

def vina_features(protname, ligname, i, return_dict):	
	'''generate the Vina features'''
	if not os.path.exists("temp_%s/%s/%s.pdbqt"%(uuid_str,ligname, ligname)):
		lig_prep('%s.mol2'%ligname, ligname, 'temp_%s/%s'%(uuid_str, ligname))
	
	cmd = '%s --receptor ../../%s.pdbqt --ligand %s.pdbqt --score_only --cpu 1' % (VINA_EXEC, protname, ligname)
	p = subprocess.Popen([cmd], shell=True, stdout= subprocess.PIPE, cwd="temp_%s/%s"%(uuid_str,ligname))
	contents = p.communicate()[0].decode("gbk").strip()	
	try:
		affinity = re.search(r'Affinity:(.*)\n',contents).group().split()[1].strip()
	except:
		affinity = None
	try:
		gauss1 = re.search(r'gauss 1(.*)\n',contents).group().split(':')[-1].strip()
	except:
		gauss1 = None	
	try:
		gauss2 = re.search(r'gauss 2(.*)\n',contents).group().split(':')[-1].strip()
	except:
		gauss2 = None	
	try:
		repulsion = re.search(r'repulsion(.*)\n',contents).group().split(':')[-1].strip()
	except:
		repulsion = None
	try:
		hydrophobic = re.search(r'hydrophobic(.*)\n',contents).group().split(':')[-1].strip()
	except:
		hydrophobic = None
	try:
		HB = re.search(r'Hydrogen(.*)$',contents).group().split(':')[-1].strip()
	except:
		HB = None	
	#return [protname, ligname, affinity, gauss1, gauss2, repulsion, hydrophobic, HB]
	return_dict[i] = [protname, ligname, affinity, gauss1, gauss2, repulsion, hydrophobic, HB]


def vina_features_process(args):
	'''generate the Vina features in parallel'''
	if not os.path.exists("%s.pdbqt"%args.protname):
		prot_prep(args.prot, args.protname)
	if not os.path.exists("lig_0/lig_0.mol2"):
		ligs_split(args.lig)
	
	lignames = [x for x in os.listdir('./temp_%s'%uuid_str) if os.path.isdir('./temp_%s/%s'%(uuid_str,x))]
	manger = Manager()
	return_dict = manger.dict()
	jobs = []
	pool = multiprocessing.Pool(multiprocessing.cpu_count())
	for i, ligname in enumerate(lignames):
		p = pool.apply_async(vina_features, args=(args.protname, ligname, i, return_dict))
		jobs.append(p)
	pool.close()
	pool.join()
	
	return return_dict.values()


def ecif_features(prot, protname, ligname, i, return_dict):
	'''generate the ECIF features'''
	if not os.path.exists("temp_%s/%s/%s.sdf"%(uuid_str, ligname, ligname)):
		cmd = "%s/obabel -imol2 %s.mol2 -osd > %s.sdf"%(BABEL_DIR, ligname, ligname)
		p = subprocess.Popen([cmd], shell=True, stdout= subprocess.PIPE, cwd="temp_%s/%s"%(uuid_str, ligname))
		p.wait()
				
	cmd = "python %s/features/ECIF_and_ELEM/ecif_calculate.py "%SCRIPT_DIR
	cmd += "--pdb_atom_keys_file %s//features/ECIF_and_ELEM/PDB_Atom_Keys.csv -p %s -l ./%s.sdf "%(SCRIPT_DIR, prot, ligname)
	cmd += " --prot_name %s --lig_name %s -d ecif"%(protname, ligname)
	p = subprocess.Popen([cmd], shell=True, stdout= subprocess.PIPE, cwd="temp_%s/%s"%(uuid_str, ligname))
	p.wait()	
	
	df = pd.read_csv('temp_%s/%s/ecif_ecif.csv'%(uuid_str, ligname), header=0, index_col=0)
	os.remove('temp_%s/%s/ecif_ecif.csv'%(uuid_str, ligname))
	#return df
	return_dict[i] = df
	


def ecif_features_process(args):
	if not os.path.exists("temp_%s/lig_0/lig_0.mol2"%uuid_str):
		ligs_split(args.lig)
	
	lignames = [x for x in os.listdir('./temp_%s'%uuid_str) if os.path.isdir('./temp_%s/%s'%(uuid_str,x))]
	manger = Manager()
	return_dict = manger.dict()
	jobs = []
	pool = multiprocessing.Pool(multiprocessing.cpu_count())
	for i, ligname in enumerate(lignames):
		p = pool.apply_async(ecif_features, args=(args.prot, args.protname, ligname, i, return_dict))
		jobs.append(p)
	pool.close()
	pool.join()
	
	return return_dict.values()

	
def obtain_features(args):
	if args.feats == 'vina':
		vina_feats = vina_features_process(args)		
		df_vina = pd.DataFrame(vina_feats)
		#df.columns = ['targets', 'pdb_id', 'lig_id', 'vina_affinity', 'vina_gauss_1', 'vina_gauss_2', 'vina_repulsion', 'vina_hydrophobic', 'vina_hydrogen']
		df_vina.columns = ['pdb_id', 'lig_id', 'vina_affinity', 'vina_gauss_1', 'vina_gauss_2', 'vina_repulsion', 'vina_hydrophobic', 'vina_hydrogen']
		df_vina.sort_values(by=['pdb_id','lig_id'], inplace=True)
		df_vina.reset_index(inplace=True)
		del df_vina['index']
		df.to_csv("%s_vina.csv"%args.outname)
		if args.bzip_output:
			os.system("bzip2 %s_vina.csv"%args.outname)
	elif args.feats == 'ecif':
		ecif_feats = ecif_features_process(args)
		df_ecif = pd.concat(ecif_feats,axis=0)
		df_ecif.sort_values(by=['pdb_id','lig_id'], inplace=True)
		df_ecif.reset_index(inplace=True)
		del df_ecif['index']
		df.to_csv("%s_ecif.csv"%args.outname)
		if args.bzip_output:
			os.system("bzip2 %s_ecif.csv"%args.outname)
	else:	
		vina_feats = vina_features_process(args)		
		df_vina = pd.DataFrame(vina_feats)
		#df.columns = ['targets', 'pdb_id', 'lig_id', 'vina_affinity', 'vina_gauss_1', 'vina_gauss_2', 'vina_repulsion', 'vina_hydrophobic', 'vina_hydrogen']
		df_vina.columns = ['pdb_id', 'lig_id', 'vina_affinity', 'vina_gauss_1', 'vina_gauss_2', 'vina_repulsion', 'vina_hydrophobic', 'vina_hydrogen']
		df_vina.sort_values(by=['pdb_id','lig_id'], inplace=True)
		df_vina.reset_index(inplace=True)
		del df_vina['index']
		df_vina.to_csv("%s_vina.csv"%args.outname)
		
		ecif_feats = ecif_features_process(args)
		df_ecif = pd.concat(ecif_feats,axis=0)
		df_ecif.sort_values(by=['pdb_id','lig_id'], inplace=True)
		df_ecif.reset_index(inplace=True)
		del df_ecif['index']
		df_ecif.to_csv("%s_ecif.csv"%args.outname)
		if args.bzip_output:
			os.system("bzip2 %s_vina.csv"%args.outname)
			os.system("bzip2 %s_ecif.csv"%args.outname)
	
	os.system("rm -rf %s.pdbqt"%args.protname)
	os.system("rm -rf temp_%s"%uuid_str)

def main():
	args = UserInput()
	obtain_features(args)
	
	
if __name__ == "__main__":
	main()




