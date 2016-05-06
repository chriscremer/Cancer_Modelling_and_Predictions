

#so I'll do leave-one-out then 2-fold (select hypers) and plot the predicted value of each one left out. 



#code for loading dataset
from load_eRNA_dataset import return_dataset_matrix

#to make roc curve
from make_roc import make_roc

#sklearn packages
from sklearn import cross_validation
from sklearn import preprocessing
from sklearn import linear_model

#other required packages
import numpy as np

#to plot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

set3_LG = ['BC060','BC005','BC030','BC026','BC050',
		'BC006','BC019','BC052','BC040','BC037',
		'BC032','BC039','BC041','BC029','BC024',
		'BC020','BC038','BC054','BC053','BC014',
		'BC059','BC051', 'BC063']

set3_HG = ['BC048','BC046','BC012','BC061','BC010',
			'BC013','BC035','BC027','BC036','BC008',
			'BC031','BC064','BC034','BC042','BC044',
			'BC002','BC017','BC022']


def get_labels(sample_names):
	''' Return an array of the labels of each sample'''

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

	X = dataset
	y = labels



	#############################
	#Leave one out, train, then predict that sample
	#############################

	#for storing the outputs of each sample and the target of that sample
	outputs = []
	targets = []

	#Leave one out cross validation
	cv_outer = cross_validation.LeaveOneOut(n=len(X))
	for samp_index, (train_index, test_index) in enumerate(cv_outer):

		#progress
		print 'Outer Fold ' + str(samp_index)

		X_train, X_test = X[train_index], X[test_index]
		y_train, y_test = y[train_index], y[test_index]

		#preprocess
		preprocessor = preprocessing.StandardScaler()
		preprocessor.fit(X_train)
		X_train = preprocessor.transform(X_train)
		X_test = preprocessor.transform(X_test)


		#####################################
		#select hyperparameters
		#####################################
		best_hyper = 0
		best_hyper_score = 0
		for hyper in [0.01, 0.1, 0.5, 1.0, 1.5, 2.0]:
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



		#####################################
		# use hyper and predict sample
		#####################################
		clf = linear_model.LogisticRegression(penalty='l1', C=best_hyper)
		output = clf.fit(X_train, y_train).predict_proba(X_test)
		score = output[0][1]

		outputs.append(score)
		targets.append(y_test[0])



	make_roc(targets, outputs, 'L1 Logistic Regression\nROC Curve', 'L1_log_reg_ROC.pdf')


	#finished tone
	print '\a'
