
import load_data
import numpy as np
import random
from sklearn.cross_validation import StratifiedKFold
from sklearn import preprocessing
from collections import Counter

#code for testing models and selecting hyperparameters by 2-fold
import hyper_select_and_test

#for LASSO
from sklearn import linear_model




if __name__ == "__main__":

	TURBT_data, TURBT_sample_order, RC_data, RC_sample_order = load_data.assemble_data('/data1/morrislab/ccremer/chemo_response/rawdata_main_genes_only.csv', '/data1/morrislab/ccremer/chemo_response/ClinpathDataBernCohort.xlsx', '/data1/morrislab/ccremer/chemo_response/SampleTranslation.xlsx')

	gene_names = load_data.get_gene_names('/data1/morrislab/ccremer/chemo_response/rawdata_main_genes_only.csv')


	pre_labels = [0]*len(TURBT_data)
	post_labels = [1]*len(RC_data)
	labels = pre_labels + post_labels
	labels = np.array(labels)

	dataset = np.concatenate((TURBT_data, RC_data), axis=0)
	sample_names = TURBT_sample_order + RC_sample_order

	X = dataset
	y = labels

	
	print X.shape
	X = X.T
	print X.shape
	X = X[[13616]]
	X = X.T
	

	all_scores = []
	all_sig_eRNAs = []
	for iteration in range(20):
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






			clf = linear_model.Lasso(alpha = 0.3)
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




	counts = Counter(all_sig_eRNAs)
	print counts

	#print 'weight 3836 ' + str(np.mean(weights_3836))
	#print 'weight 4046 ' + str(np.mean(weights_4046))
	#print 'weight 1562 ' + str(np.mean(weights_1562))


	print 'mean score ' + str(np.mean(all_scores))

	
	print gene_names[13616]




	#make matrix of rows = eRNAs
	#columns = samples seperated in to LG and HG 

	#dataset = dataset.T[[3836, 4046]]
	dataset = dataset.T[[13616]]
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



	'''
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	plt.pcolor(dataset, alpha=0.8, edgecolors='k', linewidths=0.5, cmap='PuBu_r')
	plt.ylim([0,len(dataset)])
	plt.savefig('test.pdf')
	'''



	from make_heatmap import make_heatmap_func
	make_heatmap_func(dataset.T, new_sample_names, 'pre_post_heatmap.pdf')
