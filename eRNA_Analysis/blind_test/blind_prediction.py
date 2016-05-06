


#add parent directory to path to get my packages there
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 


#code for loading dataset
from load_eRNA_dataset import return_dataset_matrix

#miscellaneous packages
import numpy as np
import random
from scipy.stats import mode
from collections import Counter

#sklearn methods
from sklearn import linear_model
from sklearn import cross_validation  
from sklearn import preprocessing

from make_heatmap import make_heatmap_func



def get_labels(sample_names):


	set3_LG = ['BC060','BC005','BC030','BC026','BC050',
			'BC006','BC019','BC052','BC040','BC037',
			'BC032','BC039','BC041','BC029','BC024',
			'BC020','BC038','BC054','BC053','BC014',
			'BC059','BC051', 'BC063']

	set3_HG = ['BC048','BC046','BC012','BC061','BC010',
				'BC013','BC035','BC027','BC036','BC008',
				'BC031','BC064','BC034','BC042','BC044',
				'BC002','BC017','BC022']


	labels = []
	for samp_name in sample_names:
		if samp_name in set3_LG:
			labels.append(0.0)
		else:
			labels.append(1.0)

	labels = np.array(labels)

	return labels


def accuracy(target, output):
	'''Model accuracy'''

	results = []
	for i in range(len(output)):
		if output[i] > 0.5 and target[i] > 0.5:
			results.append(1.0)
		elif output[i] < 0.5 and target[i] < 0.5:
			results.append(1.0)
		else:
			results.append(0.0)

	return np.mean(results)



if __name__ == "__main__":



	#############################
	#Assemble data
	#############################
	dataset, sample_names, gene_names = return_dataset_matrix('/data1/morrislab/ccremer/eRNA/bc_eRNA_1k_noOverlap_41Samples_09.txt')

	labels = get_labels(sample_names)
	X_train = dataset
	y_train = labels

	dataset, sample_names, gene_names = return_dataset_matrix('bc_eRNA_1k_noOverlap_8Samples_prediction.txt')

	X_test = dataset

	print X_train.shape
	print y_train.shape
	print X_test.shape


	#####################################
	#Scale the data
	#####################################
	preprocessor = preprocessing.StandardScaler()
	preprocessor.fit(X_train)
	X_train = preprocessor.transform(X_train)
	X_test = preprocessor.transform(X_test)



	#####################################
	#Select hyperparameters of model
	#####################################
	best_hyper = 0
	best_hyper_score = 0
	for hyper in [0.01, 0.05, 0.1, 0.5, 1.0, 2.0]:
		#accuracy sum over all crosses
		acc_sum = 0
		#2-fold cross validation for hyperparameter selection
		#cv_inner = cross_validation.StratifiedKFold(y_train, n_folds=2)
		#leave one out, takes longer but get better hyperparameters
		cv_inner = cross_validation.LeaveOneOut(n=len(X_train))
		for j, (train_index2, valid_index) in enumerate(cv_inner):

			X_train2, X_valid = X_train[train_index2], X_train[valid_index]
			y_train2, y_valid = y_train[train_index2], y_train[valid_index]

			clf = linear_model.LogisticRegression(penalty='l1', C=hyper)
			output = clf.fit(X_train2, y_train2).predict(X_valid)
			acc = accuracy(y_valid, output)

			acc_sum += acc

		#if this hyper is best then save it
		if acc_sum > best_hyper_score:
			best_hyper_score = acc_sum
			best_hyper = hyper
	print 'best hyper ' + str(best_hyper)



	#####################################
	#Use best hyper and predict test samples
	#####################################

	clf = linear_model.LogisticRegression(penalty='l1', C=best_hyper)
	output = clf.fit(X_train, y_train).predict_proba(X_test)

	discriminant_value = []
	for i in range(len(output)):
		discriminant_value.append(output[i][1])

	#print sample_names
	#print discriminant_value
	#print 
	for i in range(len(discriminant_value)):
		print sample_names[i] + ' ' + str(discriminant_value[i])
	print

	sig_feats = []
	weights = []
	for i in range(len(clf.coef_[0])):
		#print clf.coef_[0][i]
		if abs(clf.coef_[0][i]) > 0.0:
			sig_feats.append(i)
			weights.append(clf.coef_[0][i])
	#print sig_feats
	#print weights

	eRNA_names = []
	for index in sig_feats:
		eRNA_names.append(gene_names[index])
	#print eRNA_names
	#print len(eRNA_names)
	for i in range(len(eRNA_names)):
		print eRNA_names[i] + ' ' + str(weights[i])