
import xlrd
import csv
import pprint
import numpy as np
import random

from sklearn import preprocessing
from sklearn import linear_model
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn import tree
from sklearn.neighbors import KNeighborsClassifier
from sklearn import lda
from sklearn import qda
from sklearn.decomposition import PCA

import single_gene_median_analysis

import matplotlib.mlab as mlab
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from scipy import stats

from time import time

import NN_cc

from sklearn import cross_validation

from sklearn.metrics import roc_curve, auc

from scipy import interp

from sklearn.cross_validation import StratifiedKFold


def assemble_data(expr_file, clinical_file, translation_file):

	#################################################
	#extracting data from xlsx files into lists
	################################################
	sample_translation = []
	workbook = xlrd.open_workbook(translation_file)
	worksheet = workbook.sheet_by_name('Sheet1')
	num_rows = worksheet.nrows - 1
	curr_row = -1
	num_cells = worksheet.ncols - 1
	while curr_row < num_rows:
		curr_row += 1
		row = worksheet.row(curr_row)
		#print 'Row:', curr_row
		curr_cell = -1
		this_row = []
		while curr_cell < num_cells:
			curr_cell += 1
			# Cell Types: 0=Empty, 1=Text, 2=Number, 3=Date, 4=Boolean, 5=Error, 6=Blank
			cell_type = worksheet.cell_type(curr_row, curr_cell)
			cell_value = worksheet.cell_value(curr_row, curr_cell)
			#print '	', cell_type, ':', cell_value
			this_row.append(str(cell_value))
			#print cell_value

		sample_translation.append(this_row)

	clinpathdata = []
	workbook = xlrd.open_workbook(clinical_file)
	worksheet = workbook.sheet_by_name('cohort.data')
	num_rows = worksheet.nrows - 1
	curr_row = 4
	num_cells = worksheet.ncols - 1
	while curr_row < num_rows:
		curr_row += 1
		row = worksheet.row(curr_row)
		#print 'Row:', curr_row
		curr_cell = -1
		this_row = []
		while curr_cell < num_cells:
			curr_cell += 1
			# Cell Types: 0=Empty, 1=Text, 2=Number, 3=Date, 4=Boolean, 5=Error, 6=Blank
			cell_type = worksheet.cell_type(curr_row, curr_cell)
			cell_value = worksheet.cell_value(curr_row, curr_cell)
			#print '	', cell_type, ':', cell_value
			this_row.append(str(cell_value))
			#print cell_value

		clinpathdata.append(this_row)

	#print clinpathdata


	#################################################
	#list of the order of the samps in the raw data file
	################################################
	first_row = []
	with open(expr_file, 'rb') as csvfile:
		reader = csv.reader(csvfile, delimiter=',', quotechar='|')
		for row in reader:
			first_row = row
			break

	#################################################
	#selecting samples that have TURBT and that arent red
	#################################################
	red_samp_ids = ['BER104', 'BER096', 'BER062', 'BER027', 'BER020', 'BER011']

	ids = []
	response = []
	samps = []

	for row in sample_translation:
		if row[2] == 'TURBT':
			patient_numb = int(row[1].split('-')[1])

			for row2 in clinpathdata:
				#print int(float(row2[0]))
				#print patient_numb
				#print
				if int(float(row2[0])) == patient_numb:

					for samp in first_row:
						#dont use it if its in red
						if row[0] in red_samp_ids:
							break

						#make sure its in the data file
						if row[0] == samp:

							ids.append(patient_numb)
							samps.append(row[0])
							response.append(float(row2[14]) - 1.0)
							break


	'''
	print samps
	print ids
	print response
	print first_row
	print
	print len(samps)
	print len(ids)
	print len(response)
	'''

	response = np.array(response)

	return ids, response, samps, first_row


def get_all_expressions(expr_file, first_row, samples, ids, responses):

	columns = []
	name_order = []
	cur_column = 0
	for samp in first_row:
		if samp in samples:
			columns.append(cur_column)
			name_order.append(samp)
		cur_column += 1

	#gene are rows, samp are columns
	data = []
	with open(expr_file, 'rb') as csvfile:
		reader = csv.reader(csvfile, delimiter=',', quotechar='|')
		row_count = 0
		for row in reader:

			if row_count == 0:
				row_count += 1
				continue

			#if row[len(row)-1] == 'BCL2':
			this_gene = []
			for i in columns:
				this_gene.append(float(row[i]))
			data.append(this_gene)

			#for j in row:
			#	print row


			row_count += 1
			#if row_count > 40000:
			#	break

	#pprint.pprint(training_data)

	print 'Numb of genes ' + str(len(data))
	print 'Numb of samples ' + str(len(data[0]))


	data = np.array(data)
	#samples will be rows and genes are columns
	data = data.T

	return data, name_order


def get_gene_names(expr_file):

	gene_names = []

	with open(expr_file, 'rb') as csvfile:
		reader = csv.reader(csvfile, delimiter=',', quotechar='|')
		row_count = 0
		for row in reader:

			if row_count == 0:
				row_count += 1
				continue

			#gene_names.append((row[0], row[len(row)-1]))
			gene_names.append(row[len(row)-1])

			row_count += 1

	#print len(gene_names)
	#print gene_names[:1000]

	return gene_names


def filter_genes(mean_cutoff, std_cutoff, data, gene_names):


	takata_genes = ['SPRY1', 'OSBPL11', 'ZNF107', 'LIN7C', 'WDR90', 'SLC22A18', 'PNPO', 'CRKL']

	print 'Data shape prior to mean filter of : ' + str(mean_cutoff) + ' shape: ' + str(data.shape)

	#now genes will be rows
	data = data.T

	genes_to_keep = []
	new_gene_names = []
	for gene in range(len(data)):
		if np.mean(data[gene]) > mean_cutoff or gene_names[gene] in takata_genes:
			genes_to_keep.append(gene)
			new_gene_names.append(gene_names[gene])

	data = data[genes_to_keep]
	data = data.T
	gene_names = new_gene_names

	print 'Data shape post mean filter: ' + str(data.shape)

	#print 'Data shape prior to std filter of : ' + str(std_cutoff) + ' shape: ' + str(data.shape)

	#now genes will be rows
	data = data.T

	genes_to_keep = []
	new_gene_names = []
	for gene in range(len(data)):
		if np.std(data[gene]) > std_cutoff or gene_names[gene] in takata_genes:
			genes_to_keep.append(gene)
			new_gene_names.append(gene_names[gene])

	data = data[genes_to_keep]
	data = data.T
	gene_names = new_gene_names

	print 'Data shape post std filter: ' + str(data.shape)

	'''
	special_genes = ['TP53']
	#now genes will be rows
	data = data.T

	genes_to_keep = []
	new_gene_names = []
	for gene in range(len(data)):
		if gene_names[gene] in special_genes:
			genes_to_keep.append(gene)
			new_gene_names.append(gene_names[gene])

	data = data[genes_to_keep]
	data = data.T
	gene_names = new_gene_names

	print 'Data shape post special filter: ' + str(data.shape)
	'''


	return data, gene_names


def train_models2(training_data, responses, validation_sample, validation_response, gene_names, methods):

	list_of_model_training_accuraccies = []
	list_of_model_validation_accuraccies = []
	list_of_times = []

	for method in methods:
		start = time()
		clf = method
		clf.fit(training_data, responses)
		training_score = clf.score(training_data, responses)
		validation_score = clf.score([validation_sample], [validation_response])
		list_of_model_training_accuraccies.append(training_score)
		list_of_model_validation_accuraccies.append(validation_score)
		end = time()
		list_of_times.append(end - start)

	return list_of_model_training_accuraccies, list_of_model_validation_accuraccies, list_of_times


def leave_one_out(data, responses, sample_names, gene_names, methods):

	validation_accuracies = []
	training_accuracies = []
	time_list = []
	for valid_samp in range(len(data)):

		print "Validation sample : " + str(valid_samp+1) + '/' + str(len(data))

		#make a new list of the reamaining samples -> training samples
		training_data = np.delete(data, valid_samp, 0)
		responses3 = np.delete(responses, valid_samp, 0)
		#train models
		#modified_data, modified_gene_names, modified_validation_sample = L1_filter(training_data, responses3, gene_names, data[valid_samp])
		#list_of_model_training_accuraccies, list_of_model_validation_accuraccies = train_models(modified_data, responses3, modified_validation_sample, responses[valid_samp], modified_gene_names)
		
		#list_of_model_training_accuraccies, list_of_model_validation_accuraccies, list_of_times = train_models(training_data, responses3, data[valid_samp], responses[valid_samp], gene_names)
		

		list_of_model_training_accuraccies, list_of_model_validation_accuraccies, list_of_times = train_models2(training_data, responses3, data[valid_samp], responses[valid_samp], gene_names, methods)


		training_accuracies.append(list_of_model_training_accuraccies)
		validation_accuracies.append(list_of_model_validation_accuraccies)
		time_list.append(list_of_times)


	#take average over the 52 tries and save them
	training_accuracies = np.mean(np.array(training_accuracies), axis=0)
	validation_accuracies = np.mean(np.array(validation_accuracies), axis=0)
	time_list = np.mean(np.array(time_list), axis=0)

	return training_accuracies, validation_accuracies, time_list



def filter_genes2(mean_cutoff, std_cutoff, data, gene_names, valid_samp):

	#now genes will be rows
	data = data.T

	genes_to_keep = []
	new_gene_names = []
	for gene in range(len(data)):
		if np.mean(data[gene]) > mean_cutoff:
			genes_to_keep.append(gene)
			new_gene_names.append(gene_names[gene])

	data = data[genes_to_keep]
	data = data.T
	valid_samp = valid_samp.T
	valid_samp = valid_samp[genes_to_keep]
	valid_samp = valid_samp.T
	gene_names = new_gene_names

	#print 'Data shape post mean filter: ' + str(data.shape)

	#print 'Data shape prior to std filter of : ' + str(std_cutoff) + ' shape: ' + str(data.shape)

	#now genes will be rows
	data = data.T

	genes_to_keep = []
	new_gene_names = []
	for gene in range(len(data)):
		if np.std(data[gene]) > std_cutoff:
			genes_to_keep.append(gene)
			new_gene_names.append(gene_names[gene])

	data = data[genes_to_keep]
	data = data.T
	valid_samp = valid_samp.T
	valid_samp = valid_samp[genes_to_keep]
	valid_samp = valid_samp.T
	gene_names = new_gene_names


	print '# genes ' + str(len(data[0]))

	return data, gene_names, valid_samp

def L1_filter(data, responses, gene_names, validation_sample):

	
	#print 'NUmber of genes ' + str(len(data[0]))
	#print 'Number of gene names ' + str(len(gene_names))

	#Logistic Regression L1
	clf = linear_model.LogisticRegression(penalty='l1')
	clf.fit(data, responses)
	#print clf.score(data, responses)
	genes_to_keep = []
	new_gene_names = []
	for i in range(len(clf.coef_[0])):
		if clf.coef_[0][i] != 0.0:
			#print clf.coef_[0][i]
			genes_to_keep.append(i)
			new_gene_names.append(gene_names[i])
		#if clf.coef_[1][i] != 0.0:
		#	print 'what'

	data = data.T[genes_to_keep]

	data = data.T

	validation_sample = validation_sample.T[genes_to_keep]
	validation_sample = validation_sample.T

	gene_names = new_gene_names


	print '# genes ' + str(len(data[0]))
	
	return data, gene_names, validation_sample





if __name__ == "__main__":


	print 'Assembling data...'

	ids, responses, samps, first_row = assemble_data('rawdata_main_genes_only.csv', 'ClinpathDataBernCohort.xlsx', 'SampleTranslation.xlsx')

	data, sample_names = get_all_expressions('rawdata_main_genes_only.csv', first_row, samps, ids, responses)

	gene_names = get_gene_names('rawdata_main_genes_only.csv')

	print 'Data assembly complete.'


	X = data
	y = np.array(map(int, responses))

	gaussian_naive_bayes_auc = []
	mean_tpr_gnb = 0.0
	logistic_reg_L1_auc = []
	mean_tpr_lr1 = 0.0
	logistic_reg_L2_auc = []
	mean_tpr_lr2 = 0.0
	logistic_reg_elasticnet_auc = []
	mean_tpr_lrel = 0.0
	knn_auc = []
	mean_tpr_knn = 0.0
	random_forest_auc = []
	mean_tpr_rf = 0.0
	svc_auc = []
	mean_tpr_svc = 0.0
	lda_auc = []
	mean_tpr_lda = 0.0

	mean_fpr = np.linspace(0, 1, 100)
	

	for iteration in range(10):
		print 'Iter ' + str(iteration)

		random_shuffle = random.sample(range(len(y)), len(y))
		#print random_samp
		X = X[random_shuffle]
		y = y[random_shuffle]

		cv_outer = StratifiedKFold(y, n_folds=2)
		for i, (train_index, valid_index) in enumerate(cv_outer):
			#print 'Outer Fold ' + str(i)
			X_train, X_valid = X[train_index], X[valid_index]
			y_train, y_valid = y[train_index], y[valid_index]



			X_train, new_gene_names, X_valid = filter_genes2(0.1, 0.1, X_train, gene_names, X_valid)
			#X_train, new_gene_names, X_valid = L1_filter(X_train, y_train, gene_names, X_valid)



			#print list(y_train).count(0)
			#print list(y_train).count(1)

			#hyperparameter selection
			cv_inner = StratifiedKFold(y_train, n_folds=2)
			for j, (train_index2, valid_index2) in enumerate(cv_inner):
				#print 'Selecting Good Hyperparameters'
				half1, half2 = X_train[train_index2], X_train[valid_index2]
				y1, y2 = y_train[train_index2], y_train[valid_index2]

				#print list(y1).count(0)
				#print list(y1).count(1)
				#print list(y2).count(0)
				#print list(y2).count(1)

				best_log_reg_L1_auc = 0
				best_log_reg_L1_hyper = 0
				for hyper in [0.0001, 0.001, 0.01, 0.1, 1.0, 10]:
					auc_sum = 0
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

						clf = linear_model.LogisticRegression(penalty='l1', C=hyper)
						output = clf.fit(first_half, first_y).predict_proba(other_half)
						fpr, tpr, thresholds = roc_curve(other_y, output[:, 0])
						roc_auc = auc(fpr, tpr)
						auc_sum += roc_auc
					if auc_sum > best_log_reg_L1_auc:
						best_log_reg_L1_auc = auc_sum
						best_log_reg_L1_hyper = hyper
				#print 'Best L1 hyper ' + str(best_log_reg_L1_hyper)


				best_log_reg_L2_auc = 0
				best_log_reg_L2_hyper = 0
				for hyper in [0.0001, 0.001, 0.01, 0.1, 1.0, 10]:
					auc_sum = 0
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
						output = clf.fit(first_half, first_y).predict_proba(other_half)
						fpr, tpr, thresholds = roc_curve(other_y, output[:, 0])
						roc_auc = auc(fpr, tpr)
						auc_sum += roc_auc
					if auc_sum > best_log_reg_L2_auc:
						best_log_reg_L2_auc = auc_sum
						best_log_reg_L2_hyper = hyper
				#print 'Best L2 hyper ' + str(best_log_reg_L2_hyper)


				best_log_reg_el_auc = 0
				best_log_reg_el_hyper = (0,0)
				for l1_ratio in [0.0, 0.25, 0.5, 0.75, 1.0]:
					for alpha in [0.0001, 0.001, 0.01, 0.1, 1.0, 10]:
						auc_sum = 0
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

							clf = linear_model.SGDClassifier(penalty='elasticnet', l1_ratio=l1_ratio, alpha=alpha, loss='log', n_iter=5)
							output = clf.fit(first_half, first_y).predict_proba(other_half)
							fpr, tpr, thresholds = roc_curve(other_y, output[:, 0])
							roc_auc = auc(fpr, tpr)
							auc_sum += roc_auc
						if auc_sum > best_log_reg_el_auc:
							best_log_reg_el_auc = auc_sum
							best_log_reg_el_hyper = (l1_ratio, alpha)
				#print 'Best EL hyper ' + str(best_log_reg_el_hyper)

				best_knn_auc = 0
				best_knn_hyper = 0
				for hyper in range(1,5,2):
					auc_sum = 0
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

						clf = KNeighborsClassifier(n_neighbors=hyper)
						output = clf.fit(first_half, first_y).predict_proba(other_half)
						fpr, tpr, thresholds = roc_curve(other_y, output[:, 0])
						roc_auc = auc(fpr, tpr)
						auc_sum += roc_auc
					if auc_sum > best_knn_auc:
						best_knn_auc = auc_sum
						best_knn_hyper = hyper
				#print 'Best k hyper ' + str(best_knn_hyper)

				best_rf_auc = 0
				best_rf_hyper = 0
				for hyper in range(1,20,2):
					auc_sum = 0
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

						clf = RandomForestClassifier(n_estimators=hyper)
						output = clf.fit(first_half, first_y).predict_proba(other_half)
						fpr, tpr, thresholds = roc_curve(other_y, output[:, 0])
						roc_auc = auc(fpr, tpr)
						auc_sum += roc_auc
					if auc_sum > best_rf_auc:
						best_rf_auc = auc_sum
						best_rf_hyper = hyper
				#print 'Best k hyper ' + str(best_knn_hyper)

				best_svc_auc = 0
				best_svc_hyper = 0
				for hyper in [0.0001, 0.001, 0.01, 0.1, 1.0, 10]:
					auc_sum = 0
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

						clf = SVC(C=hyper, probability=True)
						output = clf.fit(first_half, first_y).predict_proba(other_half)
						fpr, tpr, thresholds = roc_curve(other_y, output[:, 0])
						roc_auc = auc(fpr, tpr)
						auc_sum += roc_auc
					if auc_sum > best_svc_auc:
						best_svc_auc = auc_sum
						best_svc_hyper = hyper
				#print 'Best L1 hyper ' + str(best_log_reg_L1_hyper)




				#since only 2 fold, I run on both halves already to get average
				#hyperparameter selection complete
				break


			#Test Models with their best hyperparameters

			clf = GaussianNB()
			output = clf.fit(X_train, y_train).predict_proba(X_valid)
			fpr, tpr, thresholds = roc_curve(y_valid, output[:, 0])
			roc_auc = auc(fpr, tpr)	
			gaussian_naive_bayes_auc.append(roc_auc)
			mean_tpr_gnb += interp(mean_fpr, fpr, tpr)
			mean_tpr_gnb[0] = 0.0

			clf = linear_model.LogisticRegression(penalty='l1', C=best_log_reg_L1_hyper)
			output = clf.fit(X_train, y_train).predict_proba(X_valid)
			fpr, tpr, thresholds = roc_curve(y_valid, output[:, 0])
			roc_auc = auc(fpr, tpr)	
			logistic_reg_L1_auc.append(roc_auc)
			mean_tpr_lr1 += interp(mean_fpr, fpr, tpr)
			mean_tpr_lr1[0] = 0.0

			clf = linear_model.LogisticRegression(penalty='l2', C=best_log_reg_L2_hyper)
			output = clf.fit(X_train, y_train).predict_proba(X_valid)
			fpr, tpr, thresholds = roc_curve(y_valid, output[:, 0])
			roc_auc = auc(fpr, tpr)	
			logistic_reg_L2_auc.append(roc_auc)
			mean_tpr_lr2 += interp(mean_fpr, fpr, tpr)
			mean_tpr_lr2[0] = 0.0

			clf = linear_model.SGDClassifier(penalty='elasticnet', l1_ratio=best_log_reg_el_hyper[0], alpha=best_log_reg_el_hyper[1], loss='log', n_iter=50)
			output = clf.fit(X_train, y_train).predict_proba(X_valid)
			fpr, tpr, thresholds = roc_curve(y_valid, output[:, 0])
			roc_auc = auc(fpr, tpr)	
			logistic_reg_elasticnet_auc.append(roc_auc)
			mean_tpr_lrel += interp(mean_fpr, fpr, tpr)
			mean_tpr_lrel[0] = 0.0

			clf = KNeighborsClassifier(n_neighbors=best_knn_hyper)
			output = clf.fit(X_train, y_train).predict_proba(X_valid)
			fpr, tpr, thresholds = roc_curve(y_valid, output[:, 0])
			roc_auc = auc(fpr, tpr)	
			knn_auc.append(roc_auc)
			mean_tpr_knn += interp(mean_fpr, fpr, tpr)
			mean_tpr_knn[0] = 0.0

			clf = RandomForestClassifier(n_estimators=best_rf_hyper)
			output = clf.fit(X_train, y_train).predict_proba(X_valid)
			fpr, tpr, thresholds = roc_curve(y_valid, output[:, 0])
			roc_auc = auc(fpr, tpr)	
			random_forest_auc.append(roc_auc)
			mean_tpr_rf += interp(mean_fpr, fpr, tpr)
			mean_tpr_rf[0] = 0.0

			clf = SVC(C=best_svc_hyper, probability=True)
			output = clf.fit(X_train, y_train).predict_proba(X_valid)
			fpr, tpr, thresholds = roc_curve(y_valid, output[:, 0])
			roc_auc = auc(fpr, tpr)	
			svc_auc.append(roc_auc)
			mean_tpr_svc += interp(mean_fpr, fpr, tpr)
			mean_tpr_svc[0] = 0.0

			clf = lda.LDA()
			output = clf.fit(X_train, y_train).predict_proba(X_valid)
			fpr, tpr, thresholds = roc_curve(y_valid, output[:, 0])
			roc_auc = auc(fpr, tpr)	
			lda_auc.append(roc_auc)
			mean_tpr_lda += interp(mean_fpr, fpr, tpr)
			mean_tpr_lda[0] = 0.0

	print 'Gaussian Naive Bayes ROC AUC =' + str(np.mean(gaussian_naive_bayes_auc))
	print 'Logistic Regression L1 ROC AUC =' + str(np.mean(logistic_reg_L1_auc))
	print 'Logistic Regression L2 ROC AUC =' + str(np.mean(logistic_reg_L2_auc))
	print 'Logistic Regression EL ROC AUC =' + str(np.mean(logistic_reg_elasticnet_auc))
	print 'KNN ROC AUC =' + str(np.mean(knn_auc))
	print 'Random Forest ROC AUC =' + str(np.mean(random_forest_auc))
	print 'SVC ROC AUC =' + str(np.mean(svc_auc))
	print 'LDA ROC AUC =' + str(np.mean(lda_auc))
	#print logistic_reg_L1_auc
	#print logistic_reg_L2_auc
	#print logistic_reg_elasticnet_auc




	mean_tpr_gnb /= len(gaussian_naive_bayes_auc)
	mean_tpr_gnb[-1] = 1.0
	mean_auc_gnb = auc(mean_fpr, mean_tpr_gnb)
	plt.plot(mean_fpr, mean_tpr_gnb, lw=1, label='Gaussian Naive Bayes (area = %0.2f)' % (mean_auc_gnb))

	mean_tpr_lr1 /= len(gaussian_naive_bayes_auc)
	mean_tpr_lr1[-1] = 1.0
	mean_auc_lr1 = auc(mean_fpr, mean_tpr_lr1)
	plt.plot(mean_fpr, mean_tpr_lr1, lw=1, label='Logistic Regression L1 (area = %0.2f)' % (mean_auc_lr1))

	mean_tpr_lr2 /= len(gaussian_naive_bayes_auc)
	mean_tpr_lr2[-1] = 1.0
	mean_auc_lr2 = auc(mean_fpr, mean_tpr_lr2)
	plt.plot(mean_fpr, mean_tpr_lr2, lw=1, label='Logistic Regression L2 (area = %0.2f)' % (mean_auc_lr2))

	mean_tpr_lrel /= len(gaussian_naive_bayes_auc)
	mean_tpr_lrel[-1] = 1.0
	mean_auc_lrel = auc(mean_fpr, mean_tpr_lrel)
	plt.plot(mean_fpr, mean_tpr_lrel, lw=1, label='Logistic Regression Elastic Net (area = %0.2f)' % (mean_auc_lrel))

	mean_tpr_knn /= len(gaussian_naive_bayes_auc)
	mean_tpr_knn[-1] = 1.0
	mean_auc_knn = auc(mean_fpr, mean_tpr_knn)
	plt.plot(mean_fpr, mean_tpr_knn, lw=1, label='k-Nearest Neighbors (area = %0.2f)' % (mean_auc_knn))

	mean_tpr_rf /= len(gaussian_naive_bayes_auc)
	mean_tpr_rf[-1] = 1.0
	mean_auc_rf = auc(mean_fpr, mean_tpr_rf)
	plt.plot(mean_fpr, mean_tpr_rf, lw=1, label='Random Forest (area = %0.2f)' % (mean_auc_rf))

	mean_tpr_svc /= len(gaussian_naive_bayes_auc)
	mean_tpr_svc[-1] = 1.0
	mean_auc_svc = auc(mean_fpr, mean_tpr_svc)
	plt.plot(mean_fpr, mean_tpr_svc, lw=1, label='Suport Vector Classifier RBF kernel (area = %0.2f)' % (mean_auc_svc))

	mean_tpr_lda /= len(gaussian_naive_bayes_auc)
	mean_tpr_lda[-1] = 1.0
	mean_auc_lda = auc(mean_fpr, mean_tpr_lda)
	plt.plot(mean_fpr, mean_tpr_lda, lw=1, label='Linear Discriminant Analysis (area = %0.2f)' % (mean_auc_lda))


	plt.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='Luck')
	#plt.plot(mean_fpr, mean_tpr_gnb, 'k--', label='Mean ROC (area = %0.2f)' % mean_auc, lw=2)

	plt.xlim([-0.05, 1.05])
	plt.ylim([-0.05, 1.05])
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	#plt.title('Receiver operating characteristic')
	plt.legend(loc="lower right", prop={'size':6})
	plt.savefig('roc_multiple_filter_10iter.pdf')
	print 'Plot complete'


	a

	#this should actually be in the nested loop, but whatever for now
	data, gene_names = filter_genes(-0.1, 0.1, data, gene_names)

	

	list_of_model_descriptions = ['LOG REG L1', 'LOG REG L2 C=0.5', 'LOG REG L2 C=1.0', 'LOG REG L2 C=2.0', 'ElasticNet', 'GNB', '3NN', '5NN', 'RandF', 'SVC', 'LDA']
	methods = [linear_model.LogisticRegression(penalty='l1'), 
				linear_model.LogisticRegression(penalty='l2', C=0.5),
				linear_model.LogisticRegression(penalty='l2', C=1.0),
				linear_model.LogisticRegression(penalty='l2', C=2.0),
				linear_model.SGDClassifier(penalty='elasticnet', l1_ratio=0.5, alpha=0.05, loss='log', n_iter=50),
				GaussianNB(),
				KNeighborsClassifier(n_neighbors=3),
				KNeighborsClassifier(n_neighbors=5),
				RandomForestClassifier(),
				SVC(),
				lda.LDA()]


	training_accuracies, validation_accuracies, time_list = leave_one_out(data, responses, sample_names, gene_names, methods)
	for model in range(len(list_of_model_descriptions)):
		print
		print 'Model: ' + list_of_model_descriptions[model] + ' | '  + str(time_list[model]) + 'sec per model'
		print 'Validation Accuracy: ' + str(validation_accuracies[model])
		print 'Training Accuracy: ' + str(training_accuracies[model])
		print


	print
	print 'Done.'
	print
