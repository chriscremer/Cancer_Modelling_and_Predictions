
#code for loading dataset
from load_eRNA_dataset import return_dataset_matrix

#code for testing models and selecting hyperparameters by 2-fold
import hyper_select_and_test
#code that prints which samples were wrong using LR2
import hyper_select_and_test_just_LR2

#miscellaneous packages
import numpy as np
import random
from scipy.stats import mode
from collections import Counter

#sklearn methods
from sklearn import linear_model

#other sklearn
from sklearn.cross_validation import StratifiedKFold
from sklearn import preprocessing
from sklearn import mixture


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


def find_metagenes(X):

	X = X.T

	#should already be preprocessed
	# preprocessor = preprocessing.StandardScaler()
	# preprocessor.fit(X)
	# X = preprocessor.transform(X)

	lowest_bic = 0
	best_numb_of_comps = 0

	
	for numb_of_comps in range(50,51,3):
		mixture_model = mixture.GMM(n_components=numb_of_comps)
		mixture_model.fit(X)
		bic_score = mixture_model.bic(X)
		if bic_score < lowest_bic or lowest_bic == 0:
			lowest_bic = bic_score
			best_numb_of_comps = numb_of_comps

	#print best_numb_of_comps
	
	#print 'numb of metas ' + str(best_numb_of_comps)
	mixture_model = mixture.GMM(n_components=best_numb_of_comps)
	mixture_model.fit(X)
	return mixture_model.predict(X), best_numb_of_comps
	

def add_metagene_features(X, metagene_list, numb_of_metas):

	new_features = []
	for samp in range(len(X)):

		samp_metas = []
		for meta in range(numb_of_metas):

			meta_sum = 0
			for i in range(len(metagene_list)):

				if meta == metagene_list[i]:
					meta_sum += X[samp][i]

			samp_metas.append(np.mean(meta_sum))

		new_features.append(samp_metas)


	new_features = np.array(new_features)
	#print new_features.shape

	#concat
	#print X.shape
	X = np.concatenate((X, new_features), axis=1)
	#print X.shape

	return X


def remove_most_features(X, numb_of_metas, list_of_ones_to_keep):

	numb_of_features = len(X[0])
	list_of_meta_indexes = range(numb_of_features-numb_of_metas, numb_of_features)
	to_keep = list_of_ones_to_keep + list_of_meta_indexes
	X = X.T
	X = X[to_keep]
	return X.T











if __name__ == "__main__":

	dataset, sample_names, gene_names = return_dataset_matrix('bc_eRNA_1k_noOverlap_41Samples_09.txt')

	labels = get_labels(sample_names)


	X = dataset
	y = labels
	sample_names2 = sample_names





	'''
	print X.shape
	X = X.T
	print X.shape
	X = X[[3836, 4046]]
	X = X.T
	'''
	
	

	#all_sig_eRNAs = []
	#weights_3836 = []
	#weights_4046 = []
	#weights_1562 = []

	all_scores = []

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

	'''
	gaussian_naive_bayes_acc = []
	lda_acc = []
	logistic_reg_L1_acc = []
	logistic_reg_L2_acc = []
	logistic_reg_elasticnet_acc = []
	knn_acc = []
	random_forest_acc = []
	svc_acc = []
	ada_acc = []
	'''

	for iteration in range(20):
		print 'Iter ' + str(iteration)

		#shuffle it up
		random_shuffle = random.sample(range(len(y)), len(y))
		X = X[random_shuffle]
		y = y[random_shuffle]
		sample_names2 = [sample_names2[i] for i in random_shuffle]

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

			#metagene_list, numb_of_metas = find_metagenes(X_train)
			#X_train = add_metagene_features(X_train, metagene_list, numb_of_metas)
			#X_test = add_metagene_features(X_test, metagene_list, numb_of_metas)

			#X_train = remove_most_features(X_train, numb_of_metas, [3836, 4046])
			#X_test = remove_most_features(X_test, numb_of_metas, [3836, 4046])

			#only use the two good eRNAs
			X_train = X_train.T[[3836, 4046]]
			X_train = X_train.T
			X_test = X_test.T[[3836, 4046]]
			X_test = X_test.T



			methods_and_accuracies = hyper_select_and_test_just_LR2.two_fold_selection(X_train, y_train, X_test, y_test, sample_names2)

			for z in range(len(methods_and_accuracies[0])):
				#if methods_and_accuracies[0][z] == 'LR2':
				#	print methods_and_accuracies[1][z]
				accuracies_dict[methods_and_accuracies[0][z]].append(methods_and_accuracies[1][z])


			'''
			clf = linear_model.Lasso(alpha = 0.05)
			clf.fit(X_train, y_train)
			predictions = clf.predict(X_test)

			score = []
			for j in range(len(predictions)):

				if predictions[j] > 0.5 and y_test[j] > 0.5:
					score.append(1.0)
				elif predictions[j] < 0.5 and y_test[j] < 0.5:
					score.append(1.0)
				else:
					score.append(0.0)

			#print np.mean(score)

			all_scores.append(np.mean(score))

			weights = clf.coef_

			sig_weights = []
			sig_eRNA = []
			for w in range(len(weights)):
				if abs(weights[w]) > 0.0: 
					sig_weights.append(weights[w])
					sig_eRNA.append(w)

			#print sig_weights
			#print sig_eRNA
			all_sig_eRNAs.extend(sig_eRNA)
			#weights_3836.append(weights[3836])
			#weights_4046.append(weights[4046])
			#weights_1562.append(weights[1562])
			'''


	for key in accuracies_dict:
		if len(accuracies_dict[key]) > 0:
			#print len(accuracies_dict[key]) 
			print key + ' ' + str(np.mean(accuracies_dict[key]))



	'''
	#mode = mode(all_sig_eRNAs)
	counts = Counter(all_sig_eRNAs)
	print counts

	#print 'weight 3836 ' + str(np.mean(weights_3836))
	#print 'weight 4046 ' + str(np.mean(weights_4046))
	#print 'weight 1562 ' + str(np.mean(weights_1562))


	print 'mean score ' + str(np.mean(all_scores))

	print gene_names[3836]
	print gene_names[4046]

	#make matrix of rows = eRNAs
	#columns = samples seperated in to LG and HG 

	dataset = dataset.T[[3836, 4046]]
	#dataset = dataset.T[[3836]]
	dataset = dataset.T

	#print dataset

	LG_HG_sorted_indexes = np.argsort(labels)

	#for i in LG_HG_sorted_indexes:
	#	print str(LG_HG_sorted_indexes[i]) + '  ' + str(labels[i])

	dataset = dataset[LG_HG_sorted_indexes]
	labels = labels[LG_HG_sorted_indexes]
	new_sample_names = []
	for i in LG_HG_sorted_indexes:
		new_sample_names.append(sample_names[i])


	#print X

	#scale all expressions

	preprocessor = preprocessing.MinMaxScaler()
	preprocessor.fit(dataset)
	dataset = preprocessor.transform(dataset)

	#print dataset
	'''

	'''
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	plt.pcolor(dataset, alpha=0.8, edgecolors='k', linewidths=0.5, cmap='PuBu_r')
	plt.ylim([0,len(dataset)])
	plt.savefig('test.pdf')
	'''


	'''
	from make_heatmap import make_heatmap_func
	make_heatmap_func(dataset.T, new_sample_names)
	'''




