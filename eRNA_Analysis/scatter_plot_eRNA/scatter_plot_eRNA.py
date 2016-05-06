



#so I'll do leave-one-out then 2-fold (select hypers) and plot the predicted value of each one left out. 



#code for loading dataset
from load_eRNA_dataset import return_dataset_matrix

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

	dataset, sample_names, gene_names = return_dataset_matrix('/data1/morrislab/ccremer/eRNA/bc_eRNA_1k_noOverlap_41Samples_09.txt')

	labels = get_labels(sample_names)


	X = dataset
	y = labels


	#onlye use the two good eRNAs
	
	print X.shape
	X = X.T
	#X = X[[3836, 4046]]
	X = X[[611, 776, 3703, 3836]]
	#X = X[[3836, 776]]
	
	X = X.T
	print X.shape
	
	
	#i forget what this does, need to lookup
	plt.vlines(1,0,1) 

	#this list is used to store the values, so that when plotting I dont overlap them
	results = []

	#Leave one out cross validation
	cv_outer = cross_validation.LeaveOneOut(n=len(X))
	for samp_index, (train_index, test_index) in enumerate(cv_outer):
		#print 'Outer Fold ' + str(samp_index)
		X_train, X_test = X[train_index], X[test_index]
		y_train, y_test = y[train_index], y[test_index]

		#preprocess
		preprocessor = preprocessing.StandardScaler()
		preprocessor.fit(X_train)
		X_train = preprocessor.transform(X_train)
		X_test = preprocessor.transform(X_test)


		#2-fold cross validation for hyperparameter selection
		cv_inner = cross_validation.StratifiedKFold(y_train, n_folds=2)
		for j, (train_index2, valid_index2) in enumerate(cv_inner):
			#print 'Selecting Good Hyperparameters'

			half1, half2 = X_train[train_index2], X_train[valid_index2]
			y1, y2 = y_train[train_index2], y_train[valid_index2]

			best_log_reg_L2_acc = 0
			best_log_reg_L2_hyper = 0
			#for hyper in [ 0.000001, 0.0001, 0.001, 0.01, 0.1, 1.0, 10]:
			#for hyper in [0.01, 0.1, 1.0, 10]:
			for hyper in [0.1, 0.5, 1.0, 1.5, 2.0]:
				acc_sum = 0
				for k in [0,1]:
					if k == 0:
						first_half = half1
						first_y = y1
						other_half = half2
						other_y = y2
					else:
						first_half = half2
						first_y = y2
						other_half = half1
						other_y= y1				

					clf = linear_model.LogisticRegression(penalty='l2', C=hyper)
					output = clf.fit(first_half, first_y).predict(other_half)
					acc = accuracy(other_y, output)
					acc_sum += acc
				if acc_sum > best_log_reg_L2_acc:
					best_log_reg_L2_acc = acc_sum
					best_log_reg_L2_hyper = hyper

			#only need one iteration since I try both halves already
			break

		#predict value of this sample
		clf = linear_model.LogisticRegression(penalty='l2', C=best_log_reg_L2_hyper)
		output = clf.fit(X_train, y_train).predict_proba(X_test)
		score = output[0][1]
		#print output

		####################################
		#time to plot
		####################################
		this_samp = sample_names[test_index]

		x_pos = 1.005
		for old_point in results:
			#if they are close together
			if np.abs(score - old_point) < 0.04:
				x_pos += 0.003

		if this_samp in set3_LG:
			plt.annotate(this_samp, xy=(x_pos, score), size=4.4, bbox=dict(boxstyle="round", fc="0.8"),  color='blue')
		elif this_samp in set3_HG:
			plt.annotate(this_samp, xy=(x_pos, score), size=4.4, bbox=dict(boxstyle="round", fc="0.8"), color='red')
		else:
			plt.annotate(this_samp, xy=(x_pos, score), size=3.9, bbox=dict(boxstyle="round", fc="0.8"), color='yellow')

		results.append(score)



	#y = numpy.ones(numpy.shape(results))   # Make all y values the same
	#plt.plot(y,results,'_',ms = 20)  # Plot a line at each location specified in a


	
	x = np.arange(1, 1.07, 0.01)
	#draws line at decision boundary (0.5)
	plt.plot(x, [0.5]*len(x), '--', c='grey')

	plt.xticks([])
	plt.xlim(1,1.07)
	plt.ylim(0,1.02)
	plt.title("Grade Predictions\n4 eRNAs")
	plt.ylabel("Score")
	plt.annotate('HG', xy=(0.95,0.515), xycoords='axes fraction', size='small', bbox=dict(boxstyle="round", fc="0.8"), color='red')
	plt.annotate('LG', xy=(0.95,0.45), xycoords='axes fraction', size='small', bbox=dict(boxstyle="round", fc="0.8"), color='blue')

	#plt.annotate('Number of Genes = ' + str(numbOfSigGenes), xy=(0.7,0.7.5), xycoords='axes fraction', size='small', bbox=dict(boxstyle="round", fc="0.8"), color='maroon')

	#plt.annotate('Testing Set = set_2', xy=(0.7,0.7), xycoords='axes fraction', size='small', bbox=dict(boxstyle="round", fc="0.8"), color='maroon')
	#plt.annotate('Training Set = set_1', xy=(0.7,0.65), xycoords='axes fraction', size='small', bbox=dict(boxstyle="round", fc="0.8"), color='maroon')



	plt.savefig('eRNA_predictions_4_NEW_ernas.pdf')


	#finished tone
	print '\a'
