

import load_data
import numpy as np
import random
from sklearn.cross_validation import StratifiedKFold
from sklearn import preprocessing

#code for testing models and selecting hyperparameters by 2-fold
import hyper_select_and_test






if __name__ == "__main__":

	TURBT_data, TURBT_gene_order, RC_data, RC_gene_order = load_data.assemble_data('/data1/morrislab/ccremer/chemo_response/rawdata_main_genes_only.csv', '/data1/morrislab/ccremer/chemo_response/ClinpathDataBernCohort.xlsx', '/data1/morrislab/ccremer/chemo_response/SampleTranslation.xlsx')

	pre_labels = [0]*len(TURBT_data)
	post_labels = [1]*len(RC_data)
	labels = pre_labels + post_labels
	labels = np.array(labels)

	dataset = np.concatenate((TURBT_data, RC_data), axis=0)

	X = dataset
	y = labels


	accuracies_dict = {}
	accuracies_dict['GNB'] = []
	accuracies_dict['LDA'] = []
	accuracies_dict['LR1'] = []
	accuracies_dict['LR2'] = []
	accuracies_dict['Elastic'] = []
	accuracies_dict['KNN'] = []
	accuracies_dict['RandForest'] = []
	accuracies_dict['SVC'] = []
	accuracies_dict['AdaBoost'] = []

	for iteration in range(4):
		print 'Iter ' + str(iteration)

		#shuffle it up
		random_shuffle = random.sample(range(len(y)), len(y))
		X = X[random_shuffle]
		y = y[random_shuffle]

		#stratified 2 fold cross validation
		cv_outer = StratifiedKFold(y, n_folds=2)
		for i, (train_index, test_index) in enumerate(cv_outer):
			#print 'Outer Fold ' + str(i)
			X_train, X_test = X[train_index], X[test_index]
			y_train, y_test = y[train_index], y[test_index]

			#preprocess
			preprocessor = preprocessing.StandardScaler()
			preprocessor.fit(X_train)
			X_train = preprocessor.transform(X_train)
			X_test = preprocessor.transform(X_test)

			methods_and_accuracies = hyper_select_and_test.two_fold_selection(X_train, y_train, X_test, y_test)

			for z in range(len(methods_and_accuracies[0])):
				accuracies_dict[methods_and_accuracies[0][z]].append(methods_and_accuracies[1][z])




	for key in accuracies_dict:
		if len(accuracies_dict[key]) > 0:
			#print len(accuracies_dict[key]) 
			print key + ' ' + str(np.mean(accuracies_dict[key]))


