#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ========================================================================
# using existing model for prediction.
# ========================================================================

import os, sys, glob, pickle, gc
import argparse
import numpy as np
from numpy.random import RandomState
import pandas as pd
import xgboost as xgb
from sklearn import preprocessing
from sklearn.metrics import roc_curve, auc, make_scorer
from sklearn import metrics 
##os.environ['CUDA_VISIBLE_DEVICES'] = "3"

def get_auc(y_test, y_pred):
    ##if the prediction_value is + , pos_label=1; else, pos_label=0
    fpr, tpr, thresholds = roc_curve(y_test, y_pred, pos_label=1)
    myauc = auc(fpr, tpr)
    return myauc


auroc_scorer = make_scorer(get_auc,greater_is_better=True, needs_proba=True)



def prepare_data(args, vtrans, scaler):
	###load the data
	#df_ref = pd.read_csv(args.ref_file, header=0, index_col=0)
	if args.features in ['ecif+vina','ecif+vina+rank']:
		df1 = pd.read_csv(args.input_file, header=0, index_col=0)
		df2 = pd.read_csv(args.input_file2, header=0, index_col=0)
		if ('vina_affinity' not in df1.columns) and ('vina_affinity' in df2.columns):
			df_temp = df1
			df1 = df2
			df2 = df_temp
			del df_temp			
		df_desc = df2['desc'].apply(lambda xs: pd.Series([int(x.strip()) for x in xs.strip('[]').split(',')]))
		if 'rank' in args.features:
			df1['rank'] = df1.lig_id.apply(lambda x: int(x.split('_')[-1])+1 if x != '1_000x' else 0)
			vina_columns = ['pdb_id', 'lig_id', 'rank','vina_affinity','vina_gauss_1','vina_gauss_2','vina_repulsion','vina_hydrophobic','vina_hydrogen']
			df1 = df1[vina_columns]	
			
		else:
			vina_columns = ['pdb_id', 'lig_id', 'vina_affinity','vina_gauss_1','vina_gauss_2','vina_repulsion','vina_hydrophobic','vina_hydrogen']
			df1 = df1[vina_columns]	
		
		df = pd.concat([df1, df_desc], axis=1)
		df['label'] = 0
		del df1, df2, df_desc
	else:
		df = pd.read_csv(args.input_file, header=0, index_col=0)
		if args.features == 'nnscore':
			#df = pd.concat([df, df_ref['label']], axis=1)
			#del df_ref
			df['label'] = 0
		elif args.features == 'nnscore-vina': 
			#df = pd.concat([df, df_ref['label']], axis=1)
			vina_columns = ['vina_affinity','vina_gauss_1','vina_gauss_2','vina_repulsion','vina_hydrophobic','vina_hydrogen']
			df.drop(vina_columns, axis=1, inplace=True)
			df['label'] = 0
			#del df_ref
		elif args.features == 'vina': 
			#df = pd.concat([df, df_ref['label']], axis=1)
			vina_columns = ['vina_affinity','vina_gauss_1','vina_gauss_2','vina_repulsion','vina_hydrophobic','vina_hydrogen']
			df = df[['pdb_id', 'lig_id'] + vina_columns + ['label']]
			#del df_ref
			df['label'] = 0
		elif args.features == 'nnscore+rank': 
			df['rank'] = df.lig_id.apply(lambda x: int(x.split('_')[-1])+1 if x != '1_000x' else 0)
			#df = pd.concat([df, df_ref['label']], axis=1)
			#del df_ref
			df['label'] = 0
		elif args.features == 'e3fp': 
			df_e3fp = df['e3fp'].apply(lambda xs: pd.Series([int(x.strip()) for x in xs.strip('[]').split(',')]))
			df = pd.concat([df[['pdb_id', 'lig_id']], df_e3fp], axis=1)
			df['label'] = 0
			#del df_ref
		else:
			df_desc = df['desc'].apply(lambda xs: pd.Series([int(x.strip()) for x in xs.strip('[]').split(',')]))
			df = pd.concat([df[['pdb_id', 'lig_id']], df_desc], axis=1)
			#del df_ref
			df['label'] = 0		
		
	
	df_test = df	
	X_test = df_test.values[:, 2:-1]
	y_test = df_test.values[:,-1]
	y_test = y_test.astype(int)
	
	##deltete the features with variance < 0.01
	X_test_vtrans = vtrans.transform(X_test)
	
	###scale the data
	X_test_standard = scaler.transform(X_test_vtrans)
	
	dtest = xgb.DMatrix(X_test_standard, label=y_test)
	
	del df, X_test, X_test_vtrans
	gc.collect()
	
	return dtest, df_test[['pdb_id','lig_id','label']]

	
def model_predict(args):
	with open(args.model,'rb') as f:
		res, vtrans, scaler = pickle.loads(f.read())
	
	dtest, df_test = prepare_data(args, vtrans, scaler)
	#if args.gpu is None:
	#	clf.set_params(tree_method='hist')
	#else:
	#	clf.set_params(tree_method='gpu_hist')
	#########it seems that the model trained on gpu can be just used on gpu rather than the cpu.
	
	y_pred = res.predict(dtest)
	#y_pred_ = clf.predict_proba(X_test)
	#y_pred = np.array([_[-1] for _ in y_pred_])
	
	df_test['pred_value'] = y_pred
	del df_test['label']
	#df_out.columns = ['pred_value']
	df_test.to_csv(args.out_file)
	df_test.sort_values(by='pred_value', ascending=False, inplace=True)
	print(df_test.iloc[:args.export_topn])
	#export_test_results(df_test.label, y_pred, args.cutoff)
	#return y_pred_proba	


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-input_file','--input_file', default='/home/shenchao/AI_pose_filter/features/crossrefined/NNscore2_refinedcrossdock.csv.bz2',
									help='Input data file (.csv), (.bz2 or .zip are supported).')
	parser.add_argument('-input_file2','--input_file2', default=None,
									help='Input data file2 (.csv), (.bz2 or .zip are supported). (default: None)')
	#parser.add_argument('-ref_file','--ref_file', default='/home/shenchao/AI_pose_filter/features/crossrefined/rmsd_statistics.csv.bz2',
	#								help='Input reference data file (.csv), (.bz2 or .zip are supported).')	
	#parser.add_argument('-core_file','--core_file', default=None)   ##'/home/shenchao/AI_pose_filter/xgb_auc2/refined_core/core_pdbid.txt'
	parser.add_argument('-model','--model', default='./final_best_model.pkl',
									help='trained model with the pickle format. (default: "./final_best_model.pkl")')	
	parser.add_argument('-f', '--features', required=True, choices = ['nnscore','nnscore-vina','vina','nnscore+rank','e3fp','elem','ecif','ecif+vina','ecif+vina+rank'], 
						help='the utilized features')
	#parser.add_argument('-t', '--type', default='core', choices=['core','casf','crosscore'],
	#								help='the type. (default: core)')  																	
	#parser.add_argument('-c', '--cutoff', type=float, default=0.5,
	#								help='Cutoff to distinguish actives from inactives. (default: 0.5.)')
	parser.add_argument('-o', '--out_file', default='xgb_out.csv')
	parser.add_argument('-en', '--export_topn', default=1, type=int, help="print the ids of the topn poses (default:1).")
	parser.add_argument('-gpu', '--gpu', default=None, type=str, help='GPU id to use.')	
	args = parser.parse_args()
	if args.gpu is not None:
		os.environ['CUDA_VISIBLE_DEVICES'] = args.gpu
	
	args = parser.parse_args()
	model_predict(args)
	
	

if __name__ == '__main__':
	main()
	
	



