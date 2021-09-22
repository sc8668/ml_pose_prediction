#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ========================================================================
# hyperparameter tuning, model training and validation.
# (the hyperparameters are tuned with the hyperopt package.)
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
from sklearn.model_selection import cross_val_score, cross_validate
from hyperopt import fmin, tpe, hp, STATUS_OK, Trials
from sklearn.model_selection import StratifiedKFold
from sklearn.feature_selection import VarianceThreshold
##os.environ['CUDA_VISIBLE_DEVICES'] = "3"

def get_auc(y_test, y_pred):
    ##if the prediction_value is + , pos_label=1; else, pos_label=0
    fpr, tpr, thresholds = roc_curve(y_test, y_pred, pos_label=1)
    myauc = auc(fpr, tpr)
    return myauc


auroc_scorer = make_scorer(get_auc,greater_is_better=True, needs_proba=True)

class BuildingModel:
	def __init__(self, args, X_train, y_train, X_test=None, y_test=None):
	#def __init__(self, args, dtrain, dtest=None):
		self.skf = StratifiedKFold(n_splits=5)
		self.X_train = X_train
		self.y_train = y_train
		self.X_test = X_test
		self.y_test = y_test		
		#self.dtrain = dtrain
		#self.dtest = dtest
		self.args = args
		self.space = {"eta": hp.loguniform("eta", np.log(0.005), np.log(0.5)),
					#"n_estimators": hp.randint("n_estimators", 300),
					"gamma": hp.uniform("gamma", 0, 1.0),  
					"max_depth": hp.randint("max_depth", 15),
					"min_child_weight": hp.randint("min_child_weight", 10),
					"subsample": hp.randint("subsample", 10),
					"colsample_bytree": hp.randint("colsample_bytree", 10),
					"colsample_bylevel": hp.randint("colsample_bylevel", 10),
					"colsample_bynode": hp.randint("colsample_bynode", 10),
					"lambda": hp.loguniform("lambda", np.log(0.001), np.log(1)),
					"alpha": hp.loguniform("alpha", np.log(0.001), np.log(1)),
					}
	
	def obtain_clf(self, params):
		if self.args.gpu is None:
			clf = xgb.XGBClassifier(objective='binary:logistic',
									eval_metric='logloss',								
									silent=1, 
									seed=self.args.random_state2,
									nthread=-1,
									**params
									)			
		else:
			clf = xgb.XGBClassifier(objective='binary:logistic',
									eval_metric='logloss',
									tree_method='gpu_hist',
									silent=1, 
									seed=self.args.random_state2,
									nthread=-1,
									**params
									)
		return clf
	
	def params_tranform(self, params):
		params["max_depth"] = params["max_depth"] + 5
		#params['n_estimators'] = params['n_estimators'] * 10 + 50
		params["subsample"] = params["subsample"] * 0.05 + 0.5
		params["colsample_bytree"] = params["colsample_bytree"] * 0.05 + 0.5
		params["colsample_bylevel"] = params["colsample_bylevel"] * 0.05 + 0.5
		params["colsample_bynode"] = params["colsample_bynode"] * 0.05 + 0.5
		params["min_child_weight"] = params["min_child_weight"] + 1
		return params 	
	
	def f(self, params):
		auroc = self.hyperopt_train_test(params)
		return {'loss': -auroc, 'status': STATUS_OK}	
		
	def hyperopt_train_test(self, params):
		params = self.params_tranform(params)
		clf = self.obtain_clf(params)
		
		auroc_list = []
		for k, (train_index, test_index) in enumerate(self.skf.split(self.X_train, self.y_train)):	
			X_tr, X_te = self.X_train[train_index], self.X_train[test_index]
			y_tr, y_te = self.y_train[train_index], self.y_train[test_index]
			dtr = xgb.DMatrix(X_tr, label=y_tr)
			dval = xgb.DMatrix(X_te, label=y_te)    
			res = xgb.train(clf.get_params(), dtr, evals=[(dval,'val')], num_boost_round=5000, early_stopping_rounds=50, verbose_eval=False)
			y_pred_proba = res.predict(dval)
			auroc = get_auc(y_te, y_pred_proba)
			auroc_list.append(auroc)
				
		return np.mean(np.array(auroc_list))

	
	def best_params_save(self, best):
		'''save the best parameters'''
		#best_params = self.params_tranform(best)
		s = pickle.dumps(best)
		with open('best_params.pkl', 'wb') as f:
			f.write(s)
	
	def best_params_load(self, tuned=True):
		if tuned:
			best_params = self.hyperparams_tuning()
		else:
			if not os.path.exists('./best_params.pkl'):
				print('the file "best_params.pkl" does not exist!')
				sys.exit(1)
			else:
				with open('best_params.pkl','rb') as f:
					best_params = pickle.loads(f.read())
		return best_params
		
	def hyperparams_tuning(self):
		trials = Trials()
		best = fmin(self.f, self.space, algo=tpe.suggest, max_evals=self.args.max_evals, trials=trials)
		best_params = self.params_tranform(best)
		print('best: %s'%best_params)
		print('loss: %s'%min(trials.losses()))
		
		self.best_params_save(best_params)
		return best_params
		
	def obtain_best_cv_scores(self, tuned=True):
		'''obtain the 5-fold CV results of the training set'''
		best_params = self.best_params_load(tuned)	
		clf = self.obtain_clf(best_params)
		
		best_iteration = 0 
		auroc_list = []
		for k, (train_index, test_index) in enumerate(self.skf.split(self.X_train, self.y_train)):	
			X_tr, X_te = self.X_train[train_index], self.X_train[test_index]
			y_tr, y_te = self.y_train[train_index], self.y_train[test_index]
			dtr = xgb.DMatrix(X_tr, label=y_tr)
			dval = xgb.DMatrix(X_te, label=y_te)    
			res = xgb.train(clf.get_params(), dtr, evals=[(dval,'val')], num_boost_round=5000, early_stopping_rounds=50, verbose_eval=False)
			best_iteration = max([best_iteration, res.best_iteration])
			y_pred_proba = res.predict(dval)
			auroc = get_auc(y_te, y_pred_proba)
			auroc_list.append(auroc)
			
		return np.array(auroc_list)	
	
	def train_and_predict(self, vtrans, scaler, tuned=False, predict=True):
		'''use the best parameters to train the model, and then predict the test set'''
		best_params = self.best_params_load(tuned=True)	
		cv_scores, best_iteration = self.obtain_best_cv_scores(tuned=True)
		export_cv_results(self.args, cv_scores)
		clf = self.obtain_clf(best_params)
		dtrain = xgb.DMatrix(self.X_train, label=self.y_train)
		res = xgb.train(clf.get_params(), dtrain, num_boost_round=best_iteration)
		#clf.fit(self.X_train, self.y_train)
		#model = pickle.dumps(clf)
		model = pickle.dumps((res, vtrans, scaler))
		with open('final_best_model.pkl', 'wb') as f:
			f.write(model)
		
		if predict:
			if self.X_test is None or self.y_test is None:
				print('the X and y of the test set should be provided!')
				sys.exit(1)			
			dtest = xgb.DMatrix(self.X_test, label=self.y_test)
			y_pred_proba = res.predict(dtest)
			#y_pred_ = clf.predict_proba(self.X_test)
			#y_pred = np.array([_[-1] for _ in y_pred_])
			return y_pred_proba
		else:
			return None



def prepare_data(args, predict=True, variance_filter=True, standard_scaler=True):
	###load the data
	df_ref = pd.read_csv(args.ref_file, header=0, index_col=0)
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
		
		df = pd.concat([df1, df_desc, df_ref['label']], axis=1)
		del df_ref, df1, df2, df_desc
	else:
		df = pd.read_csv(args.input_file, header=0, index_col=0)
		if args.features == 'nnscore':
			df = pd.concat([df, df_ref['label']], axis=1)
			del df_ref
		elif args.features == 'nnscore-vina': 
			df = pd.concat([df, df_ref['label']], axis=1)
			vina_columns = ['vina_affinity','vina_gauss_1','vina_gauss_2','vina_repulsion','vina_hydrophobic','vina_hydrogen']
			df.drop(vina_columns, axis=1, inplace=True)
			del df_ref
		elif args.features == 'vina': 
			df = pd.concat([df, df_ref['label']], axis=1)
			vina_columns = ['vina_affinity','vina_gauss_1','vina_gauss_2','vina_repulsion','vina_hydrophobic','vina_hydrogen']
			df = df[['pdb_id', 'lig_id'] + vina_columns + ['label']]
			del df_ref
		elif args.features == 'nnscore+rank': 
			df['rank'] = df.lig_id.apply(lambda x: int(x.split('_')[-1])+1 if x != '1_000x' else 0)
			df = pd.concat([df, df_ref['label']], axis=1)
			del df_ref
		elif args.features == 'e3fp': 
			df_e3fp = df['e3fp'].apply(lambda xs: pd.Series([int(x.strip()) for x in xs.strip('[]').split(',')]))
			df = pd.concat([df[['pdb_id', 'lig_id']], df_e3fp, df_ref['label']], axis=1)
			del df_ref
		else:
			df_desc = df['desc'].apply(lambda xs: pd.Series([int(x.strip()) for x in xs.strip('[]').split(',')]))
			df = pd.concat([df[['pdb_id', 'lig_id']], df_desc, df_ref['label']], axis=1)
			del df_ref		
		
	if args.validation_method == 'training_only':
		X = df.values[:, 2:-1]
		y = df.values[:, -1]
		y = y.astype(int)
		del df
	else:			
		if args.validation_method == 'refined_core':
			if args.core_file is None:
				print("if conducting refined-core splitting, the core file should be provided!")
				sys.exit(1)
			else:
				core_pdb_id = [x.strip() for x in open(args.core_file, 'r').readlines()]
				df_pdbid = df.pdb_id.drop_duplicates(keep='first')
				pdb_test = df_pdbid[df.pdb_id.isin(core_pdb_id)]
				pdb_train = df_pdbid[~df_pdbid.isin(pdb_test)]
				df_train = df[df.pdb_id.isin(pdb_train)]
				df_test = df[df.pdb_id.isin(pdb_test)]
		elif args.validation_method == 'random':
			if args.random_state is None:
				print("if conducting random splitting, the random state should be provided!")
				sys.exit(1)			
			else:		
				###train-test split(4：1 random split)
				df_pdbid = df.pdb_id.drop_duplicates(keep='first')
				pdb_train = df_pdbid.sample(frac=0.8, random_state=args.random_state)
				pdb_test = df_pdbid[~df_pdbid.isin(pdb_train)]
				df_train = df[df.pdb_id.isin(pdb_train)]
				df_test = df[df.pdb_id.isin(pdb_test)]
		elif args.validation_method == '3-CCV':
			if args.pdb_prefix is None:
				print("if conducting 3-fold, the pdb_prefix should be provided!")
				sys.exit(1)
			else:
				listnames = locals()
			for k in range(1,4):		
				listnames['pdb_list'+str(k)] = [l.strip() for l in open('%s%s'%(args.pdb_prefix, k), 'r').readlines()]			
			df_test = df[df.pdb_id.isin(listnames['pdb_list'+str(args.ccv_id)])]
			df_train = df[df.pdb_id.isin(sum([listnames['pdb_list'+str(k)] for k in range(1,4) if k != args.ccv_id],[]))]
	
		X = df_train.values[:, 2:-1]
		y = df_train.values[:, -1]
		y = y.astype(int)
		X_test = df_test.values[:, 2:-1]
		y_test = df_test.values[:, -1]
		y_test = y_test.astype(int)
		del df, df_train
		gc.collect()

	
	if variance_filter:	
		##deltete the features with variance < 0.01
		vtrans = VarianceThreshold(threshold=0.01).fit(X)
		X = vtrans.transform(X)
		if args.predict and (not args.validation_method == 'training_only'):
			X_test = vtrans.transform(X_test)
	if standard_scaler:
		###scale the data
		scaler = preprocessing.StandardScaler().fit(X)
		X = scaler.transform(X)
		if args.predict and (not args.validation_method == 'training_only'):
			X_test = scaler.transform(X_test)
	
	#dtrain = xgb.DMatrix(X, label=y)
	#if args.predict and (not args.validation_method == 'training_only'):
	#	dtest = xgb.DMatrix(X_test, label=y_test)
	
	if (not args.predict) or args.validation_method == 'training_only':
		return X, y, vtrans, scaler
		#return dtrain, vtrans, scaler
	else:
		return X, y, X_test, y_test, df_test[['pdb_id','lig_id','label']], vtrans, scaler	
		#return dtrain, dtest, df_test[['pdb_id','lig_id','label']], vtrans, scaler			


def export_cv_results(results):
	#results = return_dict.values()
	#mcc = np.array(final_values[0])
	#print('train_mcc: %s ± %s' % (mcc.mean(), mcc.std()))
	print(results)
	print('train_auroc: %s ± %s' % (results.mean(), results.std()))
	


def export_test_results(y_test, y_pred, cutoff=0.5):
	myfuc = np.frompyfunc(lambda _: 1 if _ >= cutoff else 0,1,1)
	y_pred_binary = myfuc(y_pred).astype(int)
	auroc = get_auc(y_test, y_pred)
	print('test_auroc: %s'%auroc)
	accuracy = metrics.accuracy_score(y_test, y_pred_binary)
	print('test_accuracy: %s'%accuracy)
	balanced_accuracy = metrics.balanced_accuracy_score(y_test, y_pred_binary)
	print('test_balanced_accuracy: %s'%balanced_accuracy)	
	f1 = metrics.f1_score(y_test, y_pred_binary)
	print('test_f1: %s'%f1)	
	kappa = metrics.cohen_kappa_score(y_test, y_pred_binary)
	print('test_kappa: %s'%kappa)	
	mcc = metrics.matthews_corrcoef(y_test, y_pred_binary)
	print('test_mcc: %s'%mcc)
	precision = metrics.precision_score(y_test, y_pred_binary)
	print('test_precision: %s'%precision)
	recall = metrics.recall_score(y_test, y_pred_binary)
	print('test_recall: %s'%recall)



def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-input_file','--input_file', default='/home/shenchao/AI_pose_filter/features/crossrefined/NNscore2_refinedcrossdock.csv.bz2',
									help='Input data file (.csv), (.bz2 or .zip are supported).')
	parser.add_argument('-input_file2','--input_file2', default=None,
									help='Input data file2 (.csv), (.bz2 or .zip are supported). (default: None)')
	parser.add_argument('-ref_file','--ref_file', default='/home/shenchao/AI_pose_filter/features/crossrefined/rmsd_statistics.csv.bz2',
									help='Input reference data file (.csv), (.bz2 or .zip are supported).')	
	parser.add_argument('-core_file','--core_file', default=None)   ##'/home/shenchao/AI_pose_filter/xgb_auc2/refined_core/core_pdbid.txt'
	parser.add_argument('-rand', '--random_state', type=int, default=None,
									help='Random state to conduct the train-test split. (default: None)')
	parser.add_argument('-pdb_prefix','--pdb_prefix', default=None,  ##subset
									help='The prefix of the subset file name for 3-fold CCV.') 
	parser.add_argument('-ccv_id','--ccv_id', default=1,  type=int, choices = [1, 2, 3], 
									help='the id of the subset as the test set when conducting the 3-fold CCV. (default: 1)')
	
	parser.add_argument('-f', '--features', required=True, choices = ['nnscore','nnscore-vina','vina','nnscore+rank','e3fp','elem','ecif','ecif+vina','ecif+vina+rank'], 
						help='the utilized features')
	parser.add_argument('-vm', '--validation_method', required=True, choices = ['random','3-CCV','refined_core','training_only'], 
						help='the utilized features')
	parser.add_argument('-rand2', '--random_state2', type=int, default=0,
						help='the random number seed for the model classifier.')																		
	parser.add_argument('-max_evals','--max_evals', default=100, type=int,
						help='The maximun iterations of the hyperopt. (default 100)')					
	parser.add_argument('-p', '--predict', dest='predict', action='store_true',
						help='evaluate model on test set')
	parser.add_argument('-c', '--cutoff', type=float, default=0.5,
									help='Cutoff to distinguish actives from inactives. (default: 0.5.)')
	parser.add_argument('-o', '--out_file', default='xgb_out.csv')
	parser.add_argument('-gpu', '--gpu', default=None, type=str, help='GPU id to use.')	
	args = parser.parse_args()
	if args.gpu is not None:
		os.environ['CUDA_VISIBLE_DEVICES'] = args.gpu
	
	args = parser.parse_args()
	if args.predict:	
		X, y, X_test, y_test, df_test, vtrans, scaler = prepare_data(args, predict=True, variance_filter=True, standard_scaler=True)
		mymodel = BuildingModel(args, X, y, X_test, y_test)
		#cv_scores = mymodel.obtain_best_cv_scores(tuned=True)
		#export_cv_results(cv_scores)	
		y_pred = mymodel.train_and_predict(vtrans, scaler, tuned=False, predict=True)
		df_test['pred_value'] = y_pred
		df_test.to_csv(args.out_file)	
		export_test_results(df_test.label, y_pred, args.cutoff)
	else:
		X, y, vtrans, scaler = prepare_data(args, predict=False, variance_filter=True, standard_scaler=True)
		mymodel = BuildingModel(args, X, y, X_test=None, y_test=None)
		#cv_scores = mymodel.obtain_best_cv_scores(tuned=True)
		#export_cv_results(cv_scores)
		y_pred = mymodel.train_and_predict(vtrans, scaler, tuned=False, predict=False)	
	
	

if __name__ == '__main__':
	main()
	
	



