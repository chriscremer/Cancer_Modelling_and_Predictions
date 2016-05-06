
#sklearn packages
from sklearn import cross_validation
from sklearn import preprocessing
from sklearn import linear_model


#code for splitting dataset into train and test
from arbitrary_train_test_split import split_train_test

from make_heatmap import make_heatmap_func
from make_roc import make_roc

#code for loading dataset
from load_eRNA_dataset import return_dataset_matrix

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

	X_train, y_train, X_test, y_test, names_train, names_test = split_train_test(X,y, sample_names)

	print names_train
	print names_test

	#for storing the outputs of each sample and the target of that sample
	outputs = []
	targets = []


	#####################################
	#scale
	#####################################
	preprocessor = preprocessing.StandardScaler()
	preprocessor.fit(X_train)
	X_train = preprocessor.transform(X_train)
	X_test = preprocessor.transform(X_test)


	#####################################
	#select hyperparameters
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
	# use hyper and predict samples
	#####################################

	clf = linear_model.LogisticRegression(penalty='l1', C=best_hyper)
	output = clf.fit(X_train, y_train).predict_proba(X_test)
	#print output

	output2 = []
	for i in range(len(output)):
		output2.append(output[i][1])

	print output2
	print y_test
	print accuracy(y_test, output2)


	sig_feats = []
	weights = []
	for i in range(len(clf.coef_[0])):
		#print clf.coef_[0][i]
		if abs(clf.coef_[0][i]) > 0.0:
			sig_feats.append(i)
			weights.append(clf.coef_[0][i])
	print sig_feats

	
	output3 = clf.fit(X_train, y_train).predict_proba(X_train)
	output4 = []
	for i in range(len(output3)):
		output4.append(output3[i][1])
	print accuracy(y_train, output4)




	plt.figure(2)

	for i in range(len(y_train)):
		if y_train[i] == 0.0:
			plt.scatter(1, output4[i], c='blue')
		else:
			plt.scatter(4, output4[i], c='red')

	for i in range(len(y_test)):
		if y_test[i] == 0.0:
			plt.scatter(2, output2[i], c='lightblue')
		else:
			plt.scatter(3, output2[i], c='pink')

	plt.ylabel('Prediction Score')
	plt.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off')

	plt.annotate('LG', (0.3,0.03), xycoords='figure fraction')
	plt.annotate('HG', (0.7,0.03), xycoords='figure fraction')
	plt.annotate('Training', (0.2,0.07), xycoords='figure fraction')
	plt.annotate('Test', (0.4,0.07), xycoords='figure fraction')
	plt.annotate('Test', (0.61,0.07), xycoords='figure fraction')
	plt.annotate('Training', (0.78,0.07), xycoords='figure fraction')

	#plt.xlim([-0.01, 1.01])
	plt.ylim([-0.01, 1.01])
	# x = np.arange(0, 4, 0.01)
	# plt.plot(x, [0.5]*len(x), '--', c='grey')

	plt.savefig('scatter_PS.pdf')




	#make heatmap
	eRNA_names = []
	for index in sig_feats:
		eRNA_names.append(gene_names[index])
	matrix = X.T[sig_feats]
	title = 'eRNA\nLG vs HG'
	x_ticks_names = sample_names
	y_ticks_names = eRNA_names
	output_file_name = 'eRNA_heatmap_2.pdf'
	make_heatmap_func(matrix, title, x_ticks_names, y_ticks_names, output_file_name, weights, scale_the_features=True, sort_the_samples=True, labels=labels)



	#make ROC

	make_roc(y_test, output2, 'ROC_13_test_samples.pdf')











