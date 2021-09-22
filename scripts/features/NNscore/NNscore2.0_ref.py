#!/usr/bin/env python
# -*- coding: utf-8 -*-

# =============================================================================
# the demo to calculate the features from NNscore.
# =============================================================================


import sys, os, csv, glob
import subprocess
import multiprocessing
from multiprocessing import Manager
#from multiprocessing.dummy import Pool
import pandas as pd

## set the path of module
sys.path.append('/home/shenchao/AI_pose_filter/nnscore')
from NNScore2module import PDB, binana, command_line_parameters

vina_output_list = ['vina_affinity', 'vina_gauss_1', 'vina_gauss_2', 'vina_repulsion', 'vina_hydrophobic',
                    'vina_hydrogen']
ligand_receptor_atom_type_pairs_less_than_two_half_list = ['A_MN', 'OA_SA', 'HD_N', 'N_ZN', 'A_MG', 'HD_NA', 'A_CL',
                                                           'MG_OA', 'FE_HD', 'A_OA', 'NA_ZN', 'A_N', 'C_OA', 'F_HD',
                                                           'C_HD', 'NA_SA', 'A_ZN', 'C_NA', 'N_N', 'MN_N', 'F_N',
                                                           'FE_OA', 'HD_I', 'BR_C', 'MG_NA', 'C_ZN', 'CL_MG', 'BR_OA',
                                                           'A_FE', 'CL_OA', 'CL_N', 'NA_OA', 'F_ZN', 'HD_P', 'CL_ZN',
                                                           'C_C', 'C_CL', 'FE_N', 'HD_S', 'HD_MG', 'C_F', 'A_NA',
                                                           'BR_HD', 'HD_OA', 'HD_MN', 'A_SA', 'A_F', 'HD_SA', 'A_C',
                                                           'A_A', 'F_SA', 'C_N', 'HD_ZN', 'OA_OA', 'N_SA', 'CL_FE',
                                                           'C_MN', 'CL_HD', 'OA_ZN', 'MN_OA', 'C_MG', 'F_OA', 'CD_OA',
                                                           'S_ZN', 'N_OA', 'C_SA', 'N_NA', 'A_HD', 'HD_HD', 'SA_ZN']
ligand_receptor_atom_type_pairs_less_than_four_list = ['I_N', 'OA_SA', 'FE_NA', 'HD_NA', 'A_CL', 'MG_SA', 'A_CU',
                                                       'P_SA', 'C_NA', 'MN_NA', 'F_N', 'HD_N', 'HD_I', 'CL_MG', 'HD_S',
                                                       'CL_MN', 'F_OA', 'HD_OA', 'F_HD', 'A_SA', 'A_BR', 'BR_HD',
                                                       'SA_SA', 'A_MN', 'N_ZN', 'A_MG', 'I_OA', 'C_C', 'N_S', 'N_N',
                                                       'FE_N', 'NA_SA', 'BR_N', 'MN_N', 'A_P', 'BR_C', 'A_FE', 'MN_P',
                                                       'CL_OA', 'CU_HD', 'MN_S', 'A_S', 'FE_OA', 'NA_ZN', 'P_ZN', 'A_F',
                                                       'A_C', 'A_A', 'A_N', 'HD_MN', 'A_I', 'N_SA', 'C_OA', 'MG_P',
                                                       'BR_SA', 'CU_N', 'MN_OA', 'MG_N', 'HD_HD', 'C_FE', 'CL_NA',
                                                       'MG_OA', 'A_OA', 'CL_ZN', 'BR_OA', 'HD_ZN', 'HD_P', 'OA_P',
                                                       'OA_S', 'N_P', 'A_NA', 'CL_FE', 'HD_SA', 'C_MN', 'CL_HD', 'C_MG',
                                                       'FE_HD', 'MG_S', 'NA_S', 'NA_P', 'FE_SA', 'P_S', 'C_HD', 'A_ZN',
                                                       'CL_P', 'S_SA', 'CL_S', 'OA_ZN', 'N_NA', 'MN_SA', 'CL_N',
                                                       'NA_OA', 'C_ZN', 'C_CD', 'HD_MG', 'C_F', 'C_I', 'C_CL', 'C_N',
                                                       'C_P', 'C_S', 'A_HD', 'F_SA', 'MG_NA', 'OA_OA', 'CL_SA', 'S_ZN',
                                                       'N_OA', 'C_SA', 'SA_ZN']
ligand_atom_types_list = ['A', 'C', 'CL', 'I', 'N', 'P', 'S', 'BR', 'HD', 'NA', 'F', 'OA', 'SA']
ligand_receptor_atom_type_pairs_electrostatic_list = ['I_N', 'OA_SA', 'FE_NA', 'HD_NA', 'A_CL', 'MG_SA', 'P_SA', 'C_NA',
                                                      'MN_NA', 'F_N', 'HD_N', 'HD_I', 'CL_MG', 'HD_S', 'CL_MN', 'F_OA',
                                                      'HD_OA', 'F_HD', 'A_SA', 'A_BR', 'BR_HD', 'SA_SA', 'A_MN', 'N_ZN',
                                                      'A_MG', 'I_OA', 'C_C', 'N_S', 'N_N', 'FE_N', 'NA_SA', 'BR_N',
                                                      'MN_N', 'A_P', 'BR_C', 'A_FE', 'MN_P', 'CL_OA', 'CU_HD', 'MN_S',
                                                      'A_S', 'FE_OA', 'NA_ZN', 'P_ZN', 'A_F', 'A_C', 'A_A', 'A_N',
                                                      'HD_MN', 'A_I', 'N_SA', 'C_OA', 'MG_P', 'BR_SA', 'CU_N', 'MN_OA',
                                                      'MG_N', 'HD_HD', 'C_FE', 'CL_NA', 'MG_OA', 'A_OA', 'CL_ZN',
                                                      'BR_OA', 'HD_ZN', 'HD_P', 'OA_P', 'OA_S', 'N_P', 'A_NA', 'CL_FE',
                                                      'HD_SA', 'C_MN', 'CL_HD', 'C_MG', 'FE_HD', 'MG_S', 'NA_S', 'NA_P',
                                                      'FE_SA', 'P_S', 'C_HD', 'A_ZN', 'CL_P', 'S_SA', 'CL_S', 'OA_ZN',
                                                      'N_NA', 'MN_SA', 'CL_N', 'NA_OA', 'F_ZN', 'C_ZN', 'HD_MG', 'C_F',
                                                      'C_I', 'C_CL', 'C_N', 'C_P', 'C_S', 'A_HD', 'F_SA', 'MG_NA',
                                                      'OA_OA', 'CL_SA', 'S_ZN', 'N_OA', 'C_SA', 'SA_ZN']
rotateable_bonds_count_list = ['rot_bonds']
active_site_flexibility_list = ['SIDECHAIN_OTHER', 'SIDECHAIN_ALPHA', 'BACKBONE_ALPHA', 'SIDECHAIN_BETA',
                                'BACKBONE_BETA', 'BACKBONE_OTHER']

hbonds_list = ['HDONOR-LIGAND_SIDECHAIN_BETA', 'HDONOR-LIGAND_BACKBONE_OTHER', 'HDONOR-LIGAND_SIDECHAIN_ALPHA',
               'HDONOR-RECEPTOR_SIDECHAIN_OTHER', 'HDONOR-RECEPTOR_BACKBONE_ALPHA', 'HDONOR-RECEPTOR_SIDECHAIN_BETA',
               'HDONOR-RECEPTOR_SIDECHAIN_ALPHA', 'HDONOR-LIGAND_SIDECHAIN_OTHER', 'HDONOR-LIGAND_BACKBONE_BETA',
               'HDONOR-RECEPTOR_BACKBONE_BETA', 'HDONOR-RECEPTOR_BACKBONE_OTHER', 'HDONOR-LIGAND_BACKBONE_ALPHA']

hydrophobics_list = ['SIDECHAIN_OTHER', 'SIDECHAIN_ALPHA', 'BACKBONE_ALPHA', 'SIDECHAIN_BETA', 'BACKBONE_BETA','BACKBONE_OTHER']
stacking_list = ['ALPHA', 'BETA', 'OTHER']
pi_cation_list = ['LIGAND-CHARGED_BETA', 'LIGAND-CHARGED_ALPHA', 'RECEPTOR-CHARGED_BETA', 'RECEPTOR-CHARGED_OTHER','RECEPTOR-CHARGED_ALPHA', 'LIGAND-CHARGED_OTHER']
t_shaped_list = ['ALPHA', 'BETA', 'OTHER']
salt_bridges_list = ['ALPHA', 'BETA', 'OTHER']



def get_hearder_list():
    header_list = vina_output_list + ['atp2_%s' % it for it in ligand_receptor_atom_type_pairs_less_than_two_half_list] \
                  + ['atp4_%s' % it for it in ligand_receptor_atom_type_pairs_less_than_four_list] + ['lat_%s' % it for
                                                                                                      it in
                                                                                                      ligand_atom_types_list] \
                  + ['ele_%s' % it for it in
                     ligand_receptor_atom_type_pairs_electrostatic_list] + rotateable_bonds_count_list + [
                      'siteflex_%s' % it for it in active_site_flexibility_list] \
                  + ['hbond_%s' % it for it in hbonds_list] + ['hydrophobic_%s' % it for it in hydrophobics_list] + [
                      'stacking_%s' % it for it in stacking_list] \
                  + ['pi_cation_%s' % it for it in pi_cation_list] + ['t_shaped_%s' % it for it in t_shaped_list] + [
                      'salt_bridges_%s' % it for it in salt_bridges_list]
    
    return header_list


def obtain_features(rec, lig):
    cmd = "NNScore2.py -receptor %s -ligand %s" % (rec, lig)
    params_list = cmd.split()
    cmd_params = command_line_parameters(params_list)
    receptor = PDB()
    receptor.LoadPDB_from_file(rec)
    receptor.OrigFileName = rec
    d = binana(lig, receptor, cmd_params, "", "", "")
    
    final_list = d.vina_output + d.ligand_receptor_atom_type_pairs_less_than_two_half.values() + d.ligand_receptor_atom_type_pairs_less_than_four.values() \
                 + d.ligand_atom_types.values() + d.ligand_receptor_atom_type_pairs_electrostatic.values() + d.rotateable_bonds_count.values() \
                 + d.active_site_flexibility.values() + d.hbonds.values() + d.hydrophobics.values() + d.stacking.values() + d.pi_cation.values() \
                 + d.t_shaped.values() + d.salt_bridges.values()
    
    return final_list


def write_file(output_file, outline):
    buffer = open(output_file, 'w')
    buffer.write(outline)
    buffer.close()


def prot_pre(name):
    '''prepare the protein with mgltools'''
    cmdline = 'module purge &&'
    cmdline += 'module load vina &&'
    #####remove the irons in advance because prepare_receptor4.py could not recognize the irons
    cmdline += 'cat %s_p.pdb | sed \'/HETATM/\'d > %s_p2.pdb &&'%(name, name)
    cmdline += 'prepare_receptor4.py -r %s_p2.pdb -o %s_p.pdbqt -A checkhydrogens -U nphs_lps_waters_nonstdres &&' % (name, name)
    cmdline += 'rm -rf %s_p2.pdb'%name    
    p = subprocess.Popen([cmdline], shell=True, cwd="%s/%s_prot"%(name, name))
    p.wait()


def lig_pre(name, ligname):    
    '''prepare the docking poses with mgltools'''
    cmdline = 'module purge &&'
    cmdline += 'module load vina &&'
    cmdline += 'prepare_ligand4.py -l %s.mol2 -o %s_temp.pdbqt -A checkhydrogens' % (ligname, ligname)
    p = subprocess.Popen([cmdline], shell=True, cwd="%s/%s_surflex/%s"%(name, name, ligname))
    p.wait()
    
    ####if not conducting this operation, some errors may occur for some ligands.
    lines = open('%s/%s_surflex/%s/%s_temp.pdbqt' % (name, name, ligname, ligname), 'r').readlines()
    lines_new = []
    for line in lines:
        if line.startswith('ATOM'):
            lines_new.append(line[:23]+'   ' +line[26:])
        else:
            lines_new.append(line)
    final = ''.join(lines_new)
    write_file('%s/%s_surflex/%s/%s.pdbqt' % (name, name, ligname, ligname), final)
    
    p = subprocess.Popen(["rm -rf %s_temp.pdbqt"%ligname], shell=True, cwd="%s/%s_surflex/%s"%(name, name, ligname))
    p.wait()    


def lig_pre0(name):    
    '''prepare the native poses with mgltools'''    
    cmdline = 'module purge &&'
    cmdline += 'module load vina &&'
    cmdline += 'prepare_ligand4.py -l %s_l.mol2 -o %s_temp.pdbqt -A checkhydrogens' % (name, name)
    p = subprocess.Popen([cmdline], shell=True, cwd="%s/%s_prot"%(name, name))
    p.wait()

    ####if not conducting this operation, some errors may occur for some ligands.    
    lines = open('%s/%s_prot/%s_temp.pdbqt' % (name, name, name), 'r').readlines()
    lines_new = []
    for line in lines:
        if line.startswith('ATOM'):
            lines_new.append(line[:23]+'   ' +line[26:])
        else:
            lines_new.append(line)
    final = ''.join(lines_new)
    write_file('%s/%s_prot/%s_l.pdbqt' % (name, name, name), final)
    
    p = subprocess.Popen(["rm -rf %s_temp.pdbqt"%name], shell=True, cwd="%s/%s_prot"%(name, name))
    p.wait() 


def zhenghe(name, i, return_dict):
	if not os.path.exists("%s/%s_prot/%s_p.pdbqt"%(name, name, name)):
		##generate protein.pdbqt
		prot_pre(name)
	
	if not os.path.exists("%s/%s_prot/%s_l.pdbqt"%(name, name, name)):
		##generate ligand.pdbqt
		lig_pre0(name)
	
	lignames = [x for x in os.listdir('./%s/%s_surflex'%(name, name)) if (os.path.isdir('./%s/%s_surflex/%s'%(name, name, x)) and x!='1_000x')]
	features = obtain_features('%s/%s_prot/%s_p.pdbqt' % (name, name, name), '%s/%s_prot/%s_l.pdbqt' % (name, name, name))
	features.insert(0, '1_000x')
	features.insert(0, name)
	featuress = [features]	
	for ligname in lignames:
		lig_pre(name, ligname)
		features = obtain_features('%s/%s_prot/%s_p.pdbqt' % (name, name, name),
									'%s/%s_surflex/%s/%s.pdbqt' % (name, name, ligname, ligname))
		features.insert(0, ligname)
		features.insert(0, name)
		featuress.append(features)
	return_dict[i] = featuress



def main():
	names = [x for x in os.listdir('.') if os.path.isdir(x)]
	manger = Manager()
	return_dict = manger.dict()
	jobs = []
	pool = multiprocessing.Pool(30)
	for i, name in enumerate(names):
		p = pool.apply_async(zhenghe, args=(name, i, return_dict))
		jobs.append(p)
	
	pool.close()
	pool.join()

	df = pd.DataFrame(sum(return_dict.values(), []))
	header_list = get_hearder_list()
	df.columns = ['pdb_id', 'lig_id'] + header_list
	df.sort_values(by=['pdb_id','lig_id'], inplace=True)
	df.reset_index(inplace=True)
	del df['index']
	df.to_csv("NNscore2_out.csv")
	os.system("bzip2 NNscore2_out.csv")
	



if __name__ == '__main__':
    main()
