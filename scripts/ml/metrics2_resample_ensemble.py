#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ========================================================================
# the metircs for docking power(AUC, Rp, success rate)  just for ML-based models  
# ensemble strategy for for cross-docked poses
# ========================================================================
from warnings import simplefilter
simplefilter(action='ignore', category=FutureWarning)
import os, sys, csv
import argparse
import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve, auc
import multiprocessing
from multiprocessing import Manager
from sklearn.utils import resample



def calc_success_rate2(df_groupby, score_name, pos_label, rmsd_cut=2.0, topn=1):
	'''calculate the success rate'''
	total_mol = len(df_groupby)
	if pos_label == 0:
		success_mol = df_groupby.apply(lambda x: 1 if x.rmsd.loc[x.nsmallest(topn,score_name).index].min() <= rmsd_cut else 0).sum()
	else:
		success_mol = df_groupby.apply(lambda x: 1 if x.rmsd.loc[x.nlargest(topn,score_name).index].min() <= rmsd_cut else 0).sum()
	
	return success_mol/total_mol	

	

def obtain_metircs(df_test, i, return_dict, myresample=True, score_name='pred_value', pos_label=1):	
	if myresample:
		df_test = resample(df_test, random_state=i, replace=True)
	#inter_auc = get_auc(df_test.label, df_test[score_name], pos_label)
	#if pos_label ==0:
	#	inter_rank_score = df_test[[score_name, 'rmsd']].corr('spearman').iloc[0,1]
	#else:
	#	inter_rank_score = -df_test[[score_name, 'rmsd']].corr('spearman').iloc[0,1]
	
	#df_out = pd.DataFrame(df_test.pdb_id.drop_duplicates(keep='first'))	
	#df_out.index = df_out.pdb_id
	df_groupby = df_test.groupby(by='lig')
	#df_out['intra_auc'] = df_groupby.apply(lambda x: get_auc(x.label, x[score_name], pos_label))
	#df_out['intra_rank_score'] = df_groupby.apply(lambda x: x[[score_name, 'rmsd']].corr('spearman').iloc[0,1])
	#if pos_label == 1:
	#	df_out['intra_rank_score'] = -df_out['intra_rank_score']
	
	#del df_out['pdb_id']
	#df_out.to_csv('%s_intra_target_%s.csv'%(score_name, i))
	#intra_auc = df_out.intra_auc.mean()
	#intra_rank_score = df_out.intra_rank_score.mean()
	
	#success_rate_top1 = calc_success_rate(df_groupby, score_name, pos_label, 1)
	#success_rate_top3 = calc_success_rate(df_groupby, score_name, pos_label, 3)
	#success_rate_top5 = calc_success_rate(df_groupby, score_name, pos_label, 5)
	#success_rate_top20 = calc_success_rate(df_groupby, score_name, pos_label, 20)
	sr_20_top1 = calc_success_rate2(df_groupby, score_name, pos_label, rmsd_cut=2.0, topn=1)
	sr_20_top3 = calc_success_rate2(df_groupby, score_name, pos_label, rmsd_cut=2.0, topn=3)
	sr_20_top5 = calc_success_rate2(df_groupby, score_name, pos_label, rmsd_cut=2.0, topn=5)
	sr_20_top20 = calc_success_rate2(df_groupby, score_name, pos_label, rmsd_cut=2.0, topn=20)
	#sr_10_top1 = calc_success_rate2(df_groupby, score_name, pos_label, rmsd_cut=1.0, topn=1)
	#sr_10_top3 = calc_success_rate2(df_groupby, score_name, pos_label, rmsd_cut=1.0, topn=3)
	#sr_10_top5 = calc_success_rate2(df_groupby, score_name, pos_label, rmsd_cut=1.0, topn=5)
	#sr_10_top20 = calc_success_rate2(df_groupby, score_name, pos_label, rmsd_cut=1.0, topn=20)	
	#sr_05_top1 = calc_success_rate2(df_groupby, score_name, pos_label, rmsd_cut=0.5, topn=1)
	#sr_05_top3 = calc_success_rate2(df_groupby, score_name, pos_label, rmsd_cut=0.5, topn=3)
	#sr_05_top5 = calc_success_rate2(df_groupby, score_name, pos_label, rmsd_cut=0.5, topn=5)
	#sr_05_top20 = calc_success_rate2(df_groupby, score_name, pos_label, rmsd_cut=0.5, topn=20)	
	return_dict[i] = [i, 
			#inter_auc, inter_rank_score, 
			#intra_auc, intra_rank_score, 
			#success_rate_top1, success_rate_top3, success_rate_top5, success_rate_top20,
			sr_20_top1, sr_20_top3, sr_20_top5, sr_20_top20,
			#sr_10_top1, sr_10_top3, sr_10_top5, sr_10_top20,
			#sr_05_top1, sr_05_top3, sr_05_top5, sr_05_top20
			]



def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-input_file','--input_file', required=True, ##xgb_out.csv
									help='Input data file (.csv), (.bz2 or .zip are supported).')
	parser.add_argument('-ref_file','--ref_file', default='/home/shenchao/AI_pose_filter/cross-dock_refined/resample/features/features2/rmsd_statistics.csv',
									help='Input reference data file (.csv), (.bz2 or .zip are supported).')		
	parser.add_argument('-bad_file','--bad_file', default=None)  ##/home/shenchao/AI_pose_filter/cross-dock_refined/resample/CASF/casf_crossdock/bad_combinations.txt
	parser.add_argument('-o', '--out_file', default='statistics_resample_ensemble.csv',
									help='the output file. (default: statistics_resample_ensemble.csv)')   
	#parser.add_argument('-remain_crystalposes', '--remain_crystalposes', action='store_true', default=False,
	#								help='whether to remain the crystalized poses.')   
	parser.add_argument('-t', '--type', default='crosscore-redock', choices=['crosscore-all','crosscore-redock','crosscore-cross'],
									help='the type. (default: crosscore-redock)')  
	parser.add_argument('-d', '--data_name', default='nnscore',
									choices = ['nnscore','nnscore-vina','vina','nnscore+rank','e3fp','elem','ecif','ecif+vina','ecif+vina+rank'],
									help='the name of the data utilized. (default: NNscore)')
	parser.add_argument('-resample', '--resample', action='store_true', default=False,
									help='whether to conduct the resampling.')	
	parser.add_argument('-i', '--i', default=10000, type=int,
									help='The reample times.')
	args = parser.parse_args()
	i_list = [int(x) for x in np.linspace(0,100000,args.i)]	
	if not args.resample:
		i_list = [0]
	df_ref = pd.read_csv(args.ref_file, header=0, index_col=0)
	df = pd.read_csv(args.input_file, header=0, index_col=0)
	df = pd.concat([df, df_ref.loc[df.index][['rmsd']]], axis=1)
	df['lig'] = df.pdb_id.str.split('_').apply(lambda x: x[-1])
	df['prot'] = df.pdb_id.str.split('_').apply(lambda x: x[1])
	df_test = df
	
	if args.type.startswith('crosscore'):
		
		if args.bad_file is None:
			print("bad file should be provided!")
			sys.exit(1)
			
		bad_list = [l.strip() for l in open(args.bad_file,'r').readlines()]
		if args.type == 'crosscore-redock':
			df_test = df_test[df_test.pdb_id.str.split('_').apply(lambda x: True if x[1]==x[2] else False)]
		elif args.type == 'crosscore-cross':									
			df_test = df_test[df_test.pdb_id.str.split('_').apply(lambda x: False if x[1]==x[2] else True)]
			df_test= df_test[~df_test.pdb_id.isin(bad_list)]
		else:
			df_test= df_test[~df_test.pdb_id.isin(bad_list)]
		
		
	del df
	manger = Manager()
	return_dict = manger.dict()
	pool = multiprocessing.Pool(32)
	jobs = []
	for i in i_list:
		#p = multiprocessing.Process(target=parallel_repeated, args=(i, return_dict))
		p = pool.apply_async(obtain_metircs, (df_test, i, return_dict, args.resample))
		jobs.append(p)
		#p.start()
		
	pool.close()
	pool.join()	
	
	results = return_dict.values()
	final_values = list(zip(*results))
	print("resample_%s_%s"%(args.data_name, args.type))
	print("resample_sr_2.0_top1_ensemble: %s ± %s"%(np.array(final_values[1]).mean(),np.array(final_values[1]).std()))
	print("resample_sr_2.0_top3_ensemble: %s ± %s"%(np.array(final_values[2]).mean(),np.array(final_values[2]).std()))
	print("resample_sr_2.0_top5_ensemble: %s ± %s"%(np.array(final_values[3]).mean(),np.array(final_values[3]).std()))
	print("resample_sr_2.0_top20_ensemble: %s ± %s"%(np.array(final_values[4]).mean(),np.array(final_values[4]).std()))	
	
	out_df = pd.DataFrame(results)
	out_df.columns = ['id', 
			#'inter_auc', 'inter_rank_score', 'intra_auc', 'intra_rank_score', 
			#'success_rate_top1', 'success_rate_top3', 'success_rate_top5', 'success_rate_top20',
			'sr_20_top1', 'sr_20_top3', 'sr_20_top5', 'sr_20_top20',
			#'sr_10_top1', 'sr_10_top3', 'sr_10_top5', 'sr_10_top20',
			#'sr_05_top1', 'sr_05_top3', 'sr_05_top5', 'sr_05_top20'
			]
	out_df.loc['mean'] = out_df.apply(lambda x: x.mean())
	out_df.loc['std'] = out_df.apply(lambda x: x.std())
	out_df.to_csv(args.out_file)	


	
	
if __name__ == '__main__':
    main()
	











