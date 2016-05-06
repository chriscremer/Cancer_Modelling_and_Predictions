
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


from sklearn.metrics import precision_recall_fscore_support
from sklearn.metrics import classification_report


import single_gene_median_analysis

import matplotlib.mlab as mlab
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from scipy import stats

from time import time

import NN_cc

from multiprocessing import Pool

from make_simulated_data import make_simulated



def code_inventory():


	def filter_genes(mean_cutoff, std_cutoff, data, gene_names):


		takata_genes = ['SPRY1', 'OSBPL11', 'ZNF107', 'LIN7C', 'WDR90', 'SLC22A18', 'PNPO', 'CRKL']

		#print 'Data shape prior to mean filter of : ' + str(mean_cutoff) + ' shape: ' + str(data.shape)

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

		#print 'Data shape post mean filter: ' + str(data.shape)

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

		#print 'Data shape post std filter: ' + str(data.shape)

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

		validation_sample = validation_sample[genes_to_keep]

		gene_names = new_gene_names
		

		#print 'data shape ' + str(data.shape)

		#print new_gene_names


		#print list(clf.coef_[0])
		#training_score = clf.score(training_data, responses3)
		#validation_score = clf.score([validation_sample], [validation_response])
		#list_of_model_training_accuraccies.append(training_score)
		#list_of_model_validation_accuraccies.append(validation_score)



		#return data, gene_names, validation_sample
		return genes_to_keep


	def split_into_Resp_NonResp(data, response, samp_name_order):

			row_of_responders = []
			row_of_non_responders = []
			responder_samp_names = []
			non_responder_samp_names = []
			for i in range(len(response)):
				#print response[i]
				if response[i] == 0.0:
					row_of_responders.append(i)
					responder_samp_names.append(samp_name_order[i])
				else:
					row_of_non_responders.append(i)
					non_responder_samp_names.append(samp_name_order[i])

			responders_data = data[row_of_responders]
			non_responders_data = data[row_of_non_responders]

			return responders_data, non_responders_data, responder_samp_names, non_responder_samp_names


	def nested_leave_one_out(data, responses, sample_names, gene_names):


		#data = np.ma.array(data, mask=False)
		#responses = np.ma.array(responses, mask=False)

		#data.mask[3] = True
		#responses.mask[3] = True

		test_accuracy = []
		average_validation_accuracy = []
		average_training_accuracy = []
		GNB_predictions = []
		SVC_predictions = []
		actual_response = []
		list_of_model_descriptions = ['LOG REG L1', 'LOG REG L2 C=0.5', 'LOG REG L2 C=1.0', 'LOG REG L2 C=2.0', 'ElasticNet', 'GNB', '3NN', '5NN', 'RandF', 'SVC']
		#list_of_model_descriptions = ['GNB','SVC']


		for test_samp in range(len(data)):

			print "Test sample: " + str(test_samp+1) + '/' + str(len(data))

			#print 'data shape ' + str(data.shape)

			#make a new list of the reamining samples -> validation + training samples
			train_and_validation_data = np.delete(data, test_samp, 0)
			responses2 = np.delete(responses, test_samp, 0)

			#print 'train_and_validation_data shape ' + str(train_and_validation_data.shape)
			#print 'responses2 shape ' + str(responses2.shape)


			validation_accuracies = []
			training_accuracies = []
			for valid_samp in range(len(train_and_validation_data)):

				#print "V : " + str(valid_samp+1) + '/' + str(len(train_and_validation_data))

				#make a new list of the reamaining samples -> training samples
				training_data = np.delete(train_and_validation_data, valid_samp, 0)
				responses3 = np.delete(responses2, valid_samp, 0)

				#print 'training_data shape ' + str(training_data.shape)
				#print 'responses3 shape ' + str(responses3.shape)

				modified_data, modified_gene_names, modified_validation_sample = L1_filter(training_data, responses3, gene_names, train_and_validation_data[valid_samp])
				list_of_model_training_accuraccies, list_of_model_validation_accuraccies = train_models(modified_data, responses3, modified_validation_sample, responses2[valid_samp], modified_gene_names)


				#train models
				#list_of_model_training_accuraccies, list_of_model_validation_accuraccies = train_models(training_data, responses3, train_and_validation_data[valid_samp], responses2[valid_samp], gene_names)

				#print 'len of list_of_model_training_accuraccies ' + str(len(list_of_model_training_accuraccies))
				#print 'len of list_of_model_validation_accuraccies ' + str(len(list_of_model_validation_accuraccies))


				training_accuracies.append(list_of_model_training_accuraccies)
				validation_accuracies.append(list_of_model_validation_accuraccies)

				#print 'len of training_accuracies ' + str(len(training_accuracies))
				#print 'len of validation_accuracies ' + str(len(validation_accuracies))


			#take average over the 51 tries and save them
			training_accuracies = np.array(training_accuracies)
			#print 'training_accuracies shape ' + str(training_accuracies.shape)
			training_accuracies = np.mean(training_accuracies, axis=0)
			#print 'training_accuracies shape after mean ' + str(training_accuracies.shape)
			average_training_accuracy.append(training_accuracies)

			validation_accuracies = np.mean(np.array(validation_accuracies), axis=0)
			average_validation_accuracy.append(validation_accuracies)

			#now test left out test sample
			modified_data, modified_gene_names, modified_validation_sample = L1_filter(train_and_validation_data, responses2, gene_names,  data[test_samp])
			list_of_model_training_accuraccies, list_of_model_validation_accuraccies = train_models(modified_data, responses2, modified_validation_sample, responses[test_samp], modified_gene_names)


			#list_of_model_training_accuraccies, list_of_model_validation_accuraccies = train_models(train_and_validation_data, responses2, data[test_samp], responses[test_samp], gene_names)
			#list_of_model_training_accuraccies, list_of_model_validation_accuraccies, predictions = train_models_and_return_predictions(train_and_validation_data, responses2, data[test_samp], responses[test_samp], gene_names)


			#save the result
			test_accuracy.append(list_of_model_validation_accuraccies)

			#GNB_predictions.append(predictions[0])
			#SVC_predictions.append(predictions[1])
			#actual_response.append(predictions[2])

			#print training_accuracies
			#print validation_accuracies
			#print list_of_model_validation_accuraccies

		average_training_accuracy = np.array(average_training_accuracy)
		#print 'average_training_accuracy shape ' + str(average_training_accuracy.shape)
		average_training_accuracy = np.mean(average_training_accuracy, axis=0)
		#print 'average_training_accuracy shape after mean ' + str(average_training_accuracy.shape)

		average_validation_accuracy = np.mean(np.array(average_validation_accuracy), axis=0)

		test_accuracy = np.array(test_accuracy)
		#print 'test_accuracy shape ' + str(test_accuracy.shape)
		test_accuracy = np.mean(test_accuracy, axis=0)
		#print 'test_accuracy shape after mean ' + str(test_accuracy.shape)

		#for i in range(len(GNB_predictions)):
		#	print 'GNB ' + str(GNB_predictions[i]) + 'SVC ' + str(SVC_predictions[i]) + 'actual ' + str(actual_response[i])


		return list_of_model_descriptions, test_accuracy, average_validation_accuracy, average_training_accuracy


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


	def leave_one_out(data, responses, sample_names, gene_names):

		validation_accuracies = []
		training_accuracies = []
		time_list = []
		for valid_samp in range(len(data)):

			print "Validation sample : " + str(valid_samp+1) + '/' + str(len(data))

			#make a new list of the reamaining samples -> training samples
			training_data = np.delete(data, valid_samp, 0)
			responses3 = np.delete(responses, valid_samp, 0)

			#gene_pathway = fit_pathway_model2(training_data, responses3, gene_names)
			#validation_sample = pathway_model_predict2(gene_pathway, data[valid_samp])


			#train models
			#modified_data, modified_gene_names, modified_validation_sample = L1_filter(training_data, responses3, gene_names, data[valid_samp])
			#list_of_model_training_accuraccies, list_of_model_validation_accuraccies = train_models(modified_data, responses3, modified_validation_sample, responses[valid_samp], modified_gene_names)
			
			#list_of_model_training_accuraccies, list_of_model_validation_accuraccies, list_of_times = train_models(training_data, responses3, data[valid_samp], responses[valid_samp], gene_names)
			
			methods = [linear_model.LogisticRegression(penalty='l1'), 
						linear_model.LogisticRegression(penalty='l2', C=0.5),
						linear_model.LogisticRegression(penalty='l2', C=1.0),
						linear_model.LogisticRegression(penalty='l2', C=2.0),
						linear_model.SGDClassifier(penalty='elasticnet'),
						GaussianNB(),
						KNeighborsClassifier(n_neighbors=3),
						KNeighborsClassifier(n_neighbors=5),
						RandomForestClassifier(),
						SVC(),
						qda.QDA()]

			list_of_model_training_accuraccies, list_of_model_validation_accuraccies, list_of_times = train_models2(training_data, responses3, data[valid_samp], responses[valid_samp], gene_names, methods)
			#list_of_model_training_accuraccies, list_of_model_validation_accuraccies, list_of_times = train_models2(training_data, responses3, validation_sample, responses[valid_samp], gene_names, methods)


			training_accuracies.append(list_of_model_training_accuraccies)
			validation_accuracies.append(list_of_model_validation_accuraccies)
			time_list.append(list_of_times)


		#take average over the 52 tries and save them
		training_accuracies = np.mean(np.array(training_accuracies), axis=0)
		validation_accuracies = np.mean(np.array(validation_accuracies), axis=0)
		time_list = np.mean(np.array(time_list), axis=0)

		return training_accuracies, validation_accuracies, time_list


	def leave_one_out_for_pathway(data, responses, sample_names, gene_names):

		validation_accuracies = []
		training_accuracies = []
		predictions = []

		for valid_samp in range(len(data)):

			print "Validation sample : " + str(valid_samp+1) + '/' + str(len(data))

			training_data = np.delete(data, valid_samp, 0)
			responses3 = np.delete(responses, valid_samp, 0)

			sig_genes = L1_filter(training_data, responses3, gene_names, data[valid_samp])
			#print 'data shape adter L1 filter ' + str(training_data.shape)

			gene_pathway_responder, gene_pathway_non_resp = fit_pathway_model(training_data, responses3, gene_names, sample_names, sig_genes)
			prediction = pathway_model_predict(gene_pathway_responder, gene_pathway_non_resp, data[valid_samp], sig_genes)

			if prediction == responses[valid_samp]:
				validation_accuracies.append(1.0)
			else:
				validation_accuracies.append(0.0)

			predictions.append(prediction)

			#list_of_model_training_accuraccies, list_of_model_validation_accuraccies, list_of_times = train_models2(training_data, responses3, data[valid_samp], responses[valid_samp], gene_names, methods)


		for i in range(len(predictions)):

			print str(predictions[i]) + '  ' + str(responses[i])

		#take average over the 52 tries and save them
		validation_accuracies = np.mean(validation_accuracies)


		return validation_accuracies


	def fit_pathway_model(data, responses, gene_names, samp_names, sig_genes):

		responders_data, non_responders_data, responder_samp_names, non_responder_samp_names = split_into_Resp_NonResp(data, responses, samp_names)

		#responders_data = responders_data.T
		#non_responders_data = non_responders_data.T

		gene_pathway_responder = []
		gene_pathway_non_resp = []

		for gene in sig_genes:

			training_data = np.delete(responders_data.T, gene, 0)
			training_data = training_data.T
			gene_responses = responders_data.T[gene]
			clf = linear_model.Lasso()
			clf.fit(training_data, gene_responses)
			gene_pathway_responder.append(clf)

			training_data = np.delete(non_responders_data.T, gene, 0)
			training_data = training_data.T
			gene_responses = non_responders_data.T[gene]
			clf = linear_model.Lasso()
			clf.fit(training_data, gene_responses)
			gene_pathway_non_resp.append(clf)

		return gene_pathway_responder, gene_pathway_non_resp


	def pathway_model_predict(gene_pathway_responder, gene_pathway_non_resp, validation_sample, sig_genes):

		predicted_responder_sample = []
		predicted_non_resp_sample = []

		for gene in range(len(sig_genes)):

			valid_samp = np.delete(validation_sample, sig_genes[gene], 0)

			gene_value = gene_pathway_responder[gene].predict(valid_samp)
			predicted_responder_sample.append(gene_value)
			gene_value = gene_pathway_non_resp[gene].predict(valid_samp)
			predicted_non_resp_sample.append(gene_value)

		#print 'actual ' + str(list(validation_sample))
		#print 'resp ' + str(predicted_responder_sample)
		#print 'non resp ' + str(predicted_non_resp_sample)
		#print

		responder_dif = validation_sample[sig_genes] - predicted_responder_sample
		non_resp_dif = validation_sample[sig_genes] - predicted_non_resp_sample

		responder_dif = np.linalg.norm(responder_dif)
		non_resp_dif = np.linalg.norm(non_resp_dif)

		if responder_dif < non_resp_dif:
			return 0.0
		else:
			return 1.0


	def fit_pathway_model2(data, responses, gene_names):

		#responders_data, non_responders_data, responder_samp_names, non_responder_samp_names = split_into_Resp_NonResp(data, responses, samp_names)

		#responders_data = responders_data.T
		#non_responders_data = non_responders_data.T

		#gene_pathway_responder = []
		#gene_pathway_non_resp = []

		gene_pathway = []

		for gene in range(len(data[0])):

			training_data = np.delete(data.T, gene, 0)
			training_data = training_data.T
			gene_responses = data.T[gene]
			clf = linear_model.Lasso()
			clf.fit(training_data, gene_responses)
			gene_pathway.append(clf)

			'''
			training_data = np.delete(non_responders_data.T, gene, 0)
			training_data = training_data.T
			gene_responses = non_responders_data.T[gene]
			clf = linear_model.Lasso()
			clf.fit(training_data, gene_responses)
			gene_pathway_non_resp.append(clf)
			'''

		return gene_pathway 


	def pathway_model_predict2(gene_pathway, validation_sample):

		predicted_sample = []

		for gene in range(len(validation_sample)):

			valid_samp = np.delete(validation_sample, gene, 0)

			gene_value = gene_pathway[gene].predict(valid_samp)
			predicted_sample.append(gene_value)


		#print 'actual ' + str(list(validation_sample))
		#print 'resp ' + str(predicted_responder_sample)
		#print 'non resp ' + str(predicted_non_resp_sample)
		#print

		#responder_dif = validation_sample - predicted_responder_sample
		#non_resp_dif = validation_sample - predicted_non_resp_sample

		#responder_dif = np.linalg.norm(responder_dif)
		#non_resp_dif = np.linalg.norm(non_resp_dif)

		#if responder_dif < non_resp_dif:
		#	return 0.0
		#else:
		#	return 1.0
		return predicted_sample


	def single_gene_scores(data, responses, gene_names):

		scores = []

		for test_samp in range(len(data)):

			print "Test sample: " + str(test_samp+1) + '/' + str(len(data))

			train_and_validation_data = np.delete(data, test_samp, 0)
			responses2 = np.delete(responses, test_samp, 0)

			validation_accuracies = []
			training_accuracies = []

			gene_scores = [[] for x in range(len(data[0]))]
			for valid_samp in range(len(train_and_validation_data)):
				print "			Valid sample: " + str(valid_samp+1) + '/' + str(len(train_and_validation_data))

				training_data = np.delete(train_and_validation_data, valid_samp, 0)
				responses3 = np.delete(responses2, valid_samp, 0)

				for gene in range(len(training_data[0])):

					training_data_single_gene = training_data.T[gene].T
					training_data_single_gene = np.concatenate((np.array([training_data_single_gene]), np.array([np.ones(len(training_data_single_gene))])), axis=0).T
					#print training_data_single_gene.shape

					validation_samp_single_gene = [train_and_validation_data[valid_samp].T[gene]]
					validation_samp_single_gene = np.concatenate((np.array([validation_samp_single_gene]), np.array([np.ones(len(validation_samp_single_gene))])), axis=0).T
					#print validation_samp_single_gene.shape

					clf = GaussianNB()
					clf.fit(training_data_single_gene, responses3)

					#print len([train_and_validation_data[valid_samp]])
					#print len([responses2[valid_samp]])

					score = clf.score(validation_samp_single_gene, [responses2[valid_samp]])
					#prediction = clf.predict([train_and_validation_data[valid_samp].T])
					#if prediction ==  responses2[valid_samp]:
					#	score = 1.0
					#else:
					#	score = 0.0
					gene_scores[gene].append(score)



			avg_gene_scores = []
			for i in range(len(gene_scores)):
				avg_gene_scores.append(np.mean(gene_scores[i]))

			how_well_the = sorted(avg_gene_scores, reverse=True)
			print how_well_the[:40]

			top_genes = np.argsort(avg_gene_scores)
			top_genes = top_genes[::-1]
			#print top_genes[:10]


			training_data = train_and_validation_data.T[top_genes].T
			clf = GaussianNB()
			clf.fit(training_data, responses2)

			test_sample = data[test_samp].T[top_genes].T

			score = clf.score([test_sample], [responses[test_samp]])
			scores.append(score)

		print np.mean(scores)


	def neural_net_leave_one_out_validation(data, responses):

		validation_accuracies = []
		time_list = []




		for valid_samp in range(len(data)):

			start = time()

			print "Validation sample : " + str(valid_samp+1) + '/' + str(len(data))


			#make a new list of the reamaining samples -> training samples
			training_data = np.delete(training_data, valid_samp, 0)
			responses3 = np.delete(responses, valid_samp, 0)

			training_data, new_gene_names = filter_genes(0.1, 0.1, data, gene_names)


			score = neural_net(training_data, responses3, data[valid_samp], responses[valid_samp])

			end = time()
			total = (end - start)

			validation_accuracies.append(score)
			time_list.append(total)


		#take average over the 52 tries and save them
		validation_accuracies = np.mean(validation_accuracies)
		time_list = np.mean(time_list)

		return validation_accuracies, time_list



	def split_into_train_test(ids, response, samps):

		#################################################
		#selection of training samples
		#################################################
		training_samples = []
		training_ids = []
		training_response = []
		#this plus the test must add to 15 since there are only 15 responders
		half_numb_of_training_samples = 10
		numb_of_resp = 0
		numb_of_non_resp = 0
		while numb_of_resp < half_numb_of_training_samples or numb_of_non_resp < half_numb_of_training_samples:
			if numb_of_resp < half_numb_of_training_samples:
				for i in range(len(response)):
					if response[i] == 1.0:
						training_response.append(response.pop(i))
						training_ids.append(ids.pop(i))
						training_samples.append(samps.pop(i))
						numb_of_resp += 1
						break
			if numb_of_non_resp < half_numb_of_training_samples:
				for i in range(len(response)):
					if response[i] == 2.0:
						training_response.append(response.pop(i))
						training_ids.append(ids.pop(i))
						training_samples.append(samps.pop(i))
						numb_of_non_resp += 1
						break
		'''
		print test_samples
		print test_ids
		'''
		#print training_response

	 

		#################################################
		#selection of test samples
		#################################################

		#test gets all the left over, so 5 responders and 27 non
		test_samples = samps
		test_ids = ids
		test_response = response
		'''
		test_samples = []
		test_ids = []
		test_response = []
		half_numb_of_test_samples = 7
		numb_of_resp = 0
		numb_of_non_resp = 0
		while numb_of_resp < half_numb_of_test_samples or numb_of_non_resp < half_numb_of_test_samples:
			if numb_of_resp < half_numb_of_test_samples:
				for i in range(len(response)):
					if response[i] == 1.0:
						test_response.append(response.pop(i))
						test_ids.append(ids.pop(i))
						test_samples.append(samps.pop(i))
						numb_of_resp += 1
						break
			if numb_of_non_resp < half_numb_of_test_samples:
				for i in range(len(response)):
					if response[i] == 2.0:
						test_response.append(response.pop(i))
						test_ids.append(ids.pop(i))
						test_samples.append(samps.pop(i))
						numb_of_non_resp += 1
						break
		'''
		
		#print len(test_samples)
		#print len(test_ids)
		#print len(test_response)
		

		return training_samples, training_ids, training_response, test_samples, test_ids, test_response

	def get_expressions(expr_file, first_row, training_samples, training_ids, training_response, test_samples, test_ids, test_response):

		#################################################
		#Get expressions
		#################################################
		columns_of_training_data = []
		training_samp_name_order = []
		cur_column = 0
		for samp in first_row:
			if samp in training_samples:
				columns_of_training_data.append(cur_column)
				training_samp_name_order.append(samp)
			cur_column += 1

		#print columns_of_training_data

		#gene are rows, samp are columns
		training_data = []
		with open(expr_file, 'rb') as csvfile:
			reader = csv.reader(csvfile, delimiter=',', quotechar='|')
			row_count = 0
			for row in reader:

				if row_count == 0:
					row_count += 1
					continue

				#if row[len(row)-1] == 'BCL2':
				this_gene = []
				for i in columns_of_training_data:
					this_gene.append(float(row[i]))
				training_data.append(this_gene)

				#for j in row:
				#	print row


				row_count += 1
				#if row_count > 40000:
				#	break

		#pprint.pprint(training_data)
		print 'Training Data'
		print 'Numb of genes ' + str(len(training_data))
		print 'Numb of samples ' + str(len(training_data[0]))


		#################################################
		#Assemble test data
		#################################################
		columns_of_test_data = []
		test_samp_name_order = []
		cur_column = 0
		for samp in first_row:
			if samp in test_samples:
				columns_of_test_data.append(cur_column)
				test_samp_name_order.append(samp)
			cur_column += 1

		#print columns_of_test_data

		#gene are rows, samp are columns
		test_data = []
		with open(expr_file, 'rb') as csvfile:
			reader = csv.reader(csvfile, delimiter=',', quotechar='|')
			row_count = 0
			for row in reader:

				if row_count == 0:
					row_count += 1
					continue

				#if row[len(row)-1] == 'BCL2':
				this_gene = []
				for i in columns_of_test_data:
					this_gene.append(float(row[i]))
				test_data.append(this_gene)

				#for j in row:
				#	print row


				row_count += 1
				#if row_count > 40000:
					#break

		print 'Test Data'
		print 'Numb of genes ' + str(len(test_data))
		print 'Numb of samples ' + str(len(test_data[0]))


		#################################################
		#Convert lists to arrays
		#################################################
		training_data_array = np.array(training_data)
		training_data = training_data_array.T

		training_response_array = np.array(training_response)
		training_response = training_response_array.T

		test_data_array = np.array(test_data)
		test_data = test_data_array.T

		test_response_array = np.array(test_response)
		test_response = test_response_array.T

		#samps are rows. columns are expressions
		return training_data, test_data, training_response, test_response, training_samp_name_order, test_samp_name_order

	

	def other():


		#################################################
		#Remove non significant differently expressed genes
		#################################################


		#which columns of the training data are the responders
		columns_of_responders = []
		columns_of_non_responders = []
		for i in range(len(training_response)):
			if training_response[i] == 1.0:
				columns_of_responders.append(i)
			else:
				columns_of_non_responders.append(i)

		#count number of significant genes
		count = 0
		#copy list of lists so that i can remove genes 
		#training_data2 = [list(x) for x in training_data]
		#gene are the rows of the training_data, columns are samples
		#for gene in range(len(training_data)):
		max1 = 0
		min1 = 0
		all_exps = []
		for gene in range(len(training_data)-1, -1, -1):
			resp_exp_list = []
			non_resp_exp_list = []
			for column in range(len(training_data[gene])):
				if column in columns_of_responders:
					resp_exp_list.append(training_data[gene][column])
				else:
					non_resp_exp_list.append(training_data[gene][column])

				all_exps.append(training_data[gene][column])


		print 'beginning histo'
		num_bins = 50
		# the histogram of the data
		n, bins, patches = plt.hist(all_exps, num_bins, normed=1, facecolor='green', alpha=0.5)
		# add a 'best fit' line
		#y = mlab.normpdf(bins, mu, sigma)
		#plt.plot(bins, y, 'r--')
		plt.xlabel('exp')
		plt.ylabel('prop')
		plt.title(r'Histogram of all data')

		# Tweak spacing to prevent clipping of ylabel
		plt.subplots_adjust(left=0.15)
		plt.savefig('histo.pdf')

		print 'histo complete'

		'''
		two_sample_diff_var = stats.ttest_ind(resp_exp_list, non_resp_exp_list, equal_var=False)
		#print "The t-statistic is %.3f and the p-value is %.3f." % two_sample_diff_var
		if two_sample_diff_var[1] < 0.002:
			count += 1
		else:
			training_data.pop(gene)
			test_data.pop(gene)
		'''


		print 'Numb of Sig Genes = ' + str(count)
		#training_data = training_data2
		print len(training_data)
		#print training_data

		'''
		#this is just to see the difference
		for gene in range(len(training_data)-1, -1, -1):
			resp_exp_list = []
			non_resp_exp_list = []
			for column in range(len(training_data[gene])):
				if column in columns_of_responders:
					resp_exp_list.append(training_data[gene][column])
				else:
					non_resp_exp_list.append(training_data[gene][column])

			print gene
			avg_resp_exp = np.mean(resp_exp_list)
			print avg_resp_exp
			avg_non_resp_exp = np.mean(non_resp_exp_list)
			print avg_non_resp_exp
			std_resp_exp = np.std(resp_exp_list)
			print std_resp_exp
			std_non_resp_exp = np.std(non_resp_exp_list)
			print std_non_resp_exp
		'''

	def preprocess_data(training_data, test_data, training_response, test_response):

		#################################################
		#Preprocessing -> Scaling
		#################################################

		#columns get zero mean and unit variance
		#training_data = preprocessing.scale(training_data)
		#test_data = preprocessing.scale(test_data)

		'''
		scaler = preprocessing.StandardScaler()
		scaler.fit(training_data)
		scaler.transform(training_data)
		scaler.transform(test_data)
		'''

		'''
		scaler = preprocessing.StandardScaler()
		a = np.array([[1.0,2.0],[1.0,10.0],[1.0,3.0]])
		b = np.array([[4.0,20.0],[4.0,100.0],[4.0,30.0]])
		#print scaler.mean_
		print a
		print b
		scaler.fit(a)
		#print scaler['mean_']
		scaler.transform(a)
		scaler.transform(b)
		print a
		print b
		'''

		#makes the responders 0.0 and non 1.0
		min_max_scaler = preprocessing.MinMaxScaler()
		training_response = min_max_scaler.fit_transform(training_response)
		test_response = min_max_scaler.fit_transform(test_response)

		#training_data_scaled = min_max_scaler.fit_transform(training_data)
		#try using the same scaling for test and training,..., so dont fit to test
		#test_data_scaled = min_max_scaler.transform(test_data)

		'''
		print 'beginning histo'
		num_bins = 50
		# the histogram of the data
		n, bins, patches = plt.hist(training_data.T[0], num_bins, normed=1, facecolor='green', alpha=0.5)
		# add a 'best fit' line
		#y = mlab.normpdf(bins, mu, sigma)
		#plt.plot(bins, y, 'r--')
		plt.xlabel('exp')
		plt.ylabel('prop')
		plt.title(r'Histogram of all data')

		# Tweak spacing to prevent clipping of ylabel
		plt.subplots_adjust(left=0.15)
		plt.savefig('histo_processed.pdf')

		print 'histo complete'
		'''


		#print test_data_scaled

		return training_data, test_data, training_response, test_response

	def learn_and_test(training_data, test_data, training_response, test_response):

		#################################################
		#Learn
		#################################################


		training_data_scaled = training_data
		training_response_scaled = training_response
		test_data_scaled = test_data
		test_response_scaled = test_response

		ensemble_vote_training = [[] for i in training_data_scaled]
		ensemble_vote_test = [[] for i in test_data_scaled]


		#log_reg = linear_model.LogisticRegression(C=0.001)
		clf = linear_model.LogisticRegression()
		clf.fit(training_data_scaled, training_response_scaled)
		print
		print 'Log Reg Training Score: ' + str(clf.score(training_data_scaled, training_response_scaled))
		print 'Log Reg Test Score: ' + str(clf.score(test_data_scaled, test_response_scaled))
		training_predictions = clf.predict(training_data_scaled)
		test_predictions = clf.predict(test_data_scaled)
		for i in range(len(training_predictions)):
			ensemble_vote_training[i].append(training_predictions[i])
		for i in range(len(test_predictions)):
			ensemble_vote_test[i].append(test_predictions[i])



		print
		clf = SVC()
		clf.fit(training_data_scaled, training_response_scaled) 
		print 'SVC Training Score: ' + str(clf.score(training_data_scaled, training_response_scaled))
		print 'SVC Test Score: ' + str(clf.score(test_data_scaled, test_response_scaled))
		training_predictions = clf.predict(training_data_scaled)
		test_predictions = clf.predict(test_data_scaled)
		for i in range(len(training_predictions)):
			ensemble_vote_training[i].append(training_predictions[i])
		for i in range(len(test_predictions)):
			ensemble_vote_test[i].append(test_predictions[i])

		'''
		print
		clf = SVC()
		clf.fit(training_data_array_t, training_response_array_t) 
		print 'SVC Training Non-Scaled Score: ' + str(clf.score(training_data_array_t, training_response_array_t))
		print 'SVC Test Non-Scaled Score: ' + str(clf.score(test_data_array_t, test_response_array_t))

		print
		print clf.predict(training_data_scaled)
		print training_response_scaled

		print
		print clf.predict(test_data_scaled)
		print test_response_scaled
		'''

		'''
		print
		clf = RandomForestClassifier()
		clf.fit(training_data_scaled, training_response_scaled) 
		print 'Random Forest Training Score: ' + str(clf.score(training_data_scaled, training_response_scaled))
		print 'Random Forest Test Score: ' + str(clf.score(test_data_scaled, test_response_scaled))
		training_predictions = clf.predict(training_data_scaled)
		test_predictions = clf.predict(test_data_scaled)
		for i in range(len(training_predictions)):
			ensemble_vote_training[i].append(training_predictions[i])
		for i in range(len(test_predictions)):
			ensemble_vote_test[i].append(test_predictions[i])
		'''

		'''
		print
		clf = GradientBoostingClassifier()
		clf.fit(training_data_scaled, training_response_scaled) 
		print 'GBM Training Score: ' + str(clf.score(training_data_scaled, training_response_scaled))
		print 'GBM Test Score: ' + str(clf.score(test_data_scaled, test_response_scaled))
		training_predictions = clf.predict(training_data_scaled)
		test_predictions = clf.predict(test_data_scaled)
		for i in range(len(training_predictions)):
			ensemble_vote_training[i].append(training_predictions[i])
		for i in range(len(test_predictions)):
			ensemble_vote_test[i].append(test_predictions[i])
		'''

		'''
		print
		clf = tree.DecisionTreeClassifier()
		clf.fit(training_data_scaled, training_response_scaled) 
		print 'Decision Tree Training Score: ' + str(clf.score(training_data_scaled, training_response_scaled))
		print 'Decision Tree Test Score: ' + str(clf.score(test_data_scaled, test_response_scaled))
		training_predictions = clf.predict(training_data_scaled)
		test_predictions = clf.predict(test_data_scaled)
		for i in range(len(training_predictions)):
			ensemble_vote_training[i].append(training_predictions[i])
		for i in range(len(test_predictions)):
			ensemble_vote_test[i].append(test_predictions[i])
		'''


		'''
		print
		clf = GaussianNB()
		clf.fit(training_data_scaled, training_response_scaled) 
		print 'Naive Bayes Training Score: ' + str(clf.score(training_data_scaled, training_response_scaled))
		print 'Naive Bayes Test Score: ' + str(clf.score(test_data_scaled, test_response_scaled))
		training_predictions = clf.predict(training_data_scaled)
		test_predictions = clf.predict(test_data_scaled)
		for i in range(len(training_predictions)):
			ensemble_vote_training[i].append(training_predictions[i])
		for i in range(len(test_predictions)):
			ensemble_vote_test[i].append(test_predictions[i])
		'''

		print
		clf = KNeighborsClassifier(n_neighbors=2)
		clf.fit(training_data_scaled, training_response_scaled) 
		print 'KNN Training Score: ' + str(clf.score(training_data_scaled, training_response_scaled))
		print 'KNN Test Score: ' + str(clf.score(test_data_scaled, test_response_scaled))
		training_predictions = clf.predict(training_data_scaled)
		test_predictions = clf.predict(test_data_scaled)
		for i in range(len(training_predictions)):
			ensemble_vote_training[i].append(training_predictions[i])
		for i in range(len(test_predictions)):
			ensemble_vote_test[i].append(test_predictions[i])


		#print ensemble_vote_test

		ensemble_vote_training_result = []
		for x in ensemble_vote_training:
			mode1 = stats.mode(x)
			#print 'yo'
			ensemble_vote_training_result.append(mode1[0][0])

		ensemble_vote_test_result = []
		for x in ensemble_vote_test:
			#print x
			mode2 = stats.mode(x)
			#print 'hey'
			ensemble_vote_test_result.append(mode2[0][0])

		#[int(stats.mode(x)[0]) for x in ensemble_vote_training]
		#ensemble_vote_test_result = [int(stats.mode(x)[0]) for x in ensemble_vote_test]
		print
		print np.array(ensemble_vote_training_result)
		print training_response_scaled
		print
		print np.array(ensemble_vote_test_result)
		print test_response_scaled


		training_accuracy = []
		for i in range(len(ensemble_vote_training_result)):
			if ensemble_vote_training_result[i] == training_response_scaled[i]:
				training_accuracy.append(1.0)
			else:
				training_accuracy.append(0.0)
		print 'Ensemble accuracy ' + str(np.mean(training_accuracy))

		test_accuracy = []
		for i in range(len(ensemble_vote_test_result)):
			if ensemble_vote_test_result[i] == test_response_scaled[i]:
				test_accuracy.append(1.0)
			else:
				test_accuracy.append(0.0)
		print 'Ensemble accuracy ' + str(np.mean(test_accuracy))




		#print
		#print single_gene_median_analysis.median_analysis(training_data_scaled, training_response_scaled, test_data_scaled, test_response_scaled)


		'''
		from sklearn.decomposition import PCA
		pca = PCA()
		print training_data_scaled.shape
		pca.fit(training_data_scaled)

		print
		print pca.components_.shape
		print
		print pca.explained_variance_ratio_.shape
		print pca.explained_variance_ratio_
		print 

		sorted_component_1 = sorted(pca.components_[0])
		sorted_component_1_rev = sorted(pca.components_[0], reverse=True)

		print sorted_component_1[:10]
		print sorted_component_1_rev[:10]
		'''

	def takata(training_responders, training_non_responders, test_responders, test_non_responders, gene_names, training_resp_samp_names, training_non_resp_samp_names, test_resp_samp_names, test_non_resp_samp_names):

		print 'TAKATA/KATO'

		training_responders_T = training_responders.T
		training_non_responders_T = training_non_responders.T
		test_responders_T = test_responders.T
		test_non_responders_T = test_non_responders.T



		#USING THEIR GENE SET INSTEAD OF DOING WHOLE METHOD
		#we only have 8/14 genes 
		gene_set = ['SPRY1', 'OSBPL11', 'ZNF107', 'LIN7C', 'WDR90', 'SLC22A18', 'PNPO', 'CRKL']
		index_of_genes = []
		for gene in gene_set:
			index_of_genes.append(gene_names.index(gene))

		#get expressions for just those genes in each of the groups
		train_resp_gene_set = training_responders_T[index_of_genes]
		train_non_resp_gene_set = training_non_responders_T[index_of_genes]
		test_resp_gene_set = test_responders_T[index_of_genes]
		test_non_resp_gene_set = test_non_responders_T[index_of_genes]

		#train_set = np.concatenate((train_resp_gene_set, train_non_resp_gene_set), axis=1)
		#calculate mean and std of training group for each gene
		resp_mean = []
		non_resp_mean = []
		for gene in range(len(train_resp_gene_set)):
			resp_mean.append(np.mean(train_resp_gene_set[gene]))
			non_resp_mean.append(np.mean(train_non_resp_gene_set[gene]))


		#calculate PS for training samples
		training_responders_PS = []
		for sample in range(len(training_responders)):
			responder_votes = 0
			non_responder_votes = 0
			for gene in range(len(train_resp_gene_set)):
				#V
				vote_weight = abs(train_resp_gene_set[gene][sample] - ((resp_mean[gene] + non_resp_mean[gene])/2))
				if abs(resp_mean[gene] - train_resp_gene_set[gene][sample] ) < abs(non_resp_mean[gene] - train_resp_gene_set[gene][sample]):
					responder_votes += vote_weight
				else:
					non_responder_votes += vote_weight
			#PS
			PS = ((responder_votes - non_responder_votes) / (responder_votes + non_responder_votes))
			training_responders_PS.append(PS)

		training_non_responders_PS = []
		for sample in range(len(training_non_responders)):
			responder_votes = 0
			non_responder_votes = 0
			for gene in range(len(train_non_resp_gene_set)):
				#V
				vote_weight = abs(train_non_resp_gene_set[gene][sample] - ((resp_mean[gene] + non_resp_mean[gene])/2))
				if abs(resp_mean[gene] - train_non_resp_gene_set[gene][sample] ) < abs(non_resp_mean[gene] - train_resp_gene_set[gene][sample]):
					responder_votes += vote_weight
				else:
					non_responder_votes += vote_weight
			#PS
			PS = ((responder_votes - non_responder_votes) / (responder_votes + non_responder_votes))
			training_non_responders_PS.append(PS)

		correct = []
		for guess in training_responders_PS:
			if guess > 0:
				correct.append(1)
			else:
				correct.append(0)
		for guess in training_non_responders_PS:
			if guess < 0:
				correct.append(1)
			else:
				correct.append(0)

		print 'Training Accuracy ' + str(np.mean(correct))




		#same as above but now for test set

		#means should be from training not test!!!!
		#resp_mean = []
		#non_resp_mean = []
		#for gene in range(len(test_resp_gene_set)):
		#	resp_mean.append(np.mean(test_resp_gene_set[gene]))
		#	non_resp_mean.append(np.mean(test_non_resp_gene_set[gene]))


		#calculate PS for test samples
		test_responders_PS = []
		for sample in range(len(test_responders)):
			responder_votes = 0
			non_responder_votes = 0
			for gene in range(len(test_resp_gene_set)):
				#V
				vote_weight = abs(test_resp_gene_set[gene][sample] - ((resp_mean[gene] + non_resp_mean[gene])/2))
				if abs(resp_mean[gene] - test_resp_gene_set[gene][sample] ) < abs(non_resp_mean[gene] - test_resp_gene_set[gene][sample]):
					responder_votes += vote_weight
				else:
					non_responder_votes += vote_weight
			#PS
			PS = ((responder_votes - non_responder_votes) / (responder_votes + non_responder_votes))
			test_responders_PS.append(PS)

		test_non_responders_PS = []
		for sample in range(len(test_non_responders)):
			responder_votes = 0
			non_responder_votes = 0
			for gene in range(len(test_non_resp_gene_set)):
				#V
				vote_weight = abs(test_non_resp_gene_set[gene][sample] - ((resp_mean[gene] + non_resp_mean[gene])/2))
				if abs(resp_mean[gene] - test_non_resp_gene_set[gene][sample] ) < abs(non_resp_mean[gene] - test_non_resp_gene_set[gene][sample]):
					responder_votes += vote_weight
				else:
					non_responder_votes += vote_weight
			#PS
			PS = ((responder_votes - non_responder_votes) / (responder_votes + non_responder_votes))
			test_non_responders_PS.append(PS)

		correct = []
		for guess in test_responders_PS:
			if guess > 0:
				correct.append(1)
			else:
				correct.append(0)
		for guess in test_non_responders_PS:
			if guess < 0:
				correct.append(1)
			else:
				correct.append(0)

		print 'Test Accuracy ' + str(np.mean(correct))


		
		plt.figure(2)

		plt.scatter(np.ones(len(training_responders_PS)), training_responders_PS, c='blue')
		plt.scatter(np.ones(len(test_responders_PS))*2, test_responders_PS, c='lightblue')
		plt.scatter(np.ones(len(training_non_responders_PS))*4, training_non_responders_PS, c='red')	
		plt.scatter(np.ones(len(test_non_responders_PS))*3, test_non_responders_PS, c='pink')

		plt.title('Prediction Scores')
		plt.ylabel('PS')
		plt.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off')

		plt.annotate('Responders', (0.26,0.07), xycoords='figure fraction')
		plt.annotate('Non-Responders', (0.64,0.07), xycoords='figure fraction')
		plt.annotate('Training', (0.2,0.03), xycoords='figure fraction')
		plt.annotate('Test', (0.4,0.03), xycoords='figure fraction')
		plt.annotate('Test', (0.61,0.03), xycoords='figure fraction')
		plt.annotate('Training', (0.8,0.03), xycoords='figure fraction')


		plt.savefig('scatter_PS.pdf')
		


		'''
		#calculate DS for each gene
		resp_mean = []
		resp_std = []
		for gene in range(len(training_responders_T)):
			resp_mean.append(np.mean(training_responders_T[gene]))
			resp_std.append(np.std(training_responders_T[gene]))
		non_resp_mean = []
		non_resp_std = []
		for gene in range(len(training_non_responders_T)):
			non_resp_mean.append(np.mean(training_non_responders_T[gene]))
			non_resp_std.append(np.std(training_non_responders_T[gene]))
		DS = []
		for gene in range(len(training_responders_T)):
			distinguishing_score = (resp_mean[gene] - non_resp_mean[gene]) / (resp_std[gene] + non_resp_std[gene])
			DS.append(distinguishing_score)

		#random permutation test
		print 'random permutations'

		random_DS = []
		for i in range(len(DS)):
			random_DS.append([])

		numb_training_samples = len(training_responders)
		training_samples_indexes = range(numb_training_samples)
		#number of permutations
		for i in range(1000):
			print i

			#randomly select half of the training samples
			indexes_of_responders = random.sample(training_samples_indexes,numb_training_samples/2)
			indexes_of_non_responders = []
			for j in range(len(training_samples_indexes)):
				if j not in indexes_of_responders:
					indexes_of_non_responders.append(j)

			#make new sets
			random_responders = []
			for j in indexes_of_responders:
				#if it is an index for a responder
				if j < numb_training_samples/2:
					random_responders.append(training_responders[j])
				else:
					random_responders.append(training_non_responders[j-(numb_training_samples/2)])

			random_non_responders = []
			for j in indexes_of_non_responders:
				if j < numb_training_samples/2:
					random_non_responders.append(training_responders[j])
				else:
					random_non_responders.append(training_non_responders[j-(numb_training_samples/2)])

			random_responders = np.array(random_responders)
			random_non_responders = np.array(random_non_responders)

			resp_mean = []
			resp_std = []
			non_resp_mean = []
			non_resp_std = []
			for gene in range(len(random_responders.T)):
				resp_mean.append(np.mean(random_responders.T[gene]))
				resp_std.append(np.std(random_responders.T[gene]))
				non_resp_mean.append(np.mean(random_non_responders.T[gene]))
				non_resp_std.append(np.std(random_non_responders.T[gene]))
				#compute DS for each gene
				distinguishing_score = (resp_mean[gene] - non_resp_mean[gene]) / (resp_std[gene] + non_resp_std[gene])
				random_DS[gene].append(distinguishing_score)
		'''


		'''
		plt.figure(44)
		num_bins = 50
		# the histogram of the data
		n, bins, patches = plt.hist(random_DS[0], num_bins, normed=1, facecolor='green', alpha=0.5)
		# add a 'best fit' line
		#y = mlab.normpdf(bins, mu, sigma)
		#plt.plot(bins, y, 'r--')
		plt.xlabel('DS')
		plt.ylabel('Proportion')
		plt.title(r'Histogram of Distinguishing Scores')
		# Tweak spacing to prevent clipping of ylabel
		plt.subplots_adjust(left=0.15)
		plt.savefig('histo_normal_DS.pdf')
		'''

		'''
		#see if orginal DS is significant 
		rand_DS_mean = 0
		rand_DS_std = 0
		zscore = 0
		p_values = []
		up_or_down = []
		for gene in range(len(random_DS)):
			rand_DS_mean = np.mean(random_DS[gene])
			rand_DS_std = np.std(random_DS[gene])
			zscore = (DS[gene] - rand_DS_mean)/rand_DS_std
			if zscore > 0:
				p_values.append(stats.norm.sf(zscore))
				up_or_down.append(1)
			else:
				p_values.append(1 - stats.norm.sf(zscore))
				up_or_down.append(0)
		'''


		'''
		plt.figure(2)
		num_bins = 50
		# the histogram of the data
		n, bins, patches = plt.hist(p_values, num_bins, normed=1, facecolor='green', alpha=0.5)
		# add a 'best fit' line
		#y = mlab.normpdf(bins, mu, sigma)
		#plt.plot(bins, y, 'r--')
		plt.xlabel('p values')
		plt.ylabel('Proportion')
		plt.title(r'Histogram of DS p values')
		# Tweak spacing to prevent clipping of ylabel
		plt.subplots_adjust(left=0.15)
		plt.savefig('histo_pvalues.pdf')
		'''

		'''
		sort_index = np.argsort(p_values)
		sorted_pvalues = sorted(p_values)
		top_gene_names = []
		dis_gene_matrix_resp = np.array([training_responders_T[sort_index[0]]])
		dis_gene_matrix_non = np.array([training_non_responders_T[sort_index[0]]])
		for i in range(1,12):
			top_gene_names.append(gene_names[sort_index[i]])
			print gene_names[sort_index[i]]
			print sorted_pvalues[i]
			if up_or_down[sort_index[i]] == 1:
				dis_gene_matrix_resp = np.concatenate((dis_gene_matrix_resp, [training_responders_T[sort_index[i]]]), axis=0)
				dis_gene_matrix_non = np.concatenate((dis_gene_matrix_non, [training_non_responders_T[sort_index[i]]]), axis=0)
			else:
				dis_gene_matrix_resp = np.concatenate(([training_responders_T[sort_index[i]]], dis_gene_matrix_resp), axis=0)
				dis_gene_matrix_non = np.concatenate(([training_non_responders_T[sort_index[i]]], dis_gene_matrix_non), axis=0)



		print dis_gene_matrix_resp.shape
		print dis_gene_matrix_non.shape

		'''

		#matrix using only the sig genes
		#dis_gene_matrix_resp = training_responders_T[sort_index[:12]]
		#print dis_gene_matrix_resp.shape 
		#dis_gene_matrix_non = training_non_responders_T[sort_index[:12]]
		#print dis_gene_matrix_non.shape
		'''
		both_matrix = np.concatenate((dis_gene_matrix_resp, dis_gene_matrix_non), axis=1)
		print both_matrix.shape	
		'''
		#print both_matrix

		#I want the range to be centered at 0 so the heatmap is black at 0
		#max1 = both_matrix.max()
		#min1 = both_matrix.min()
		#center = (max1+min1)/2.0




		'''
		for row in range(len(both_matrix)):
			median = np.median(both_matrix[row])
			#print median
			for expr in range(len(both_matrix[row])):
				both_matrix[row][expr] = both_matrix[row][expr] - median
		'''

		#print both_matrix

		'''
		#same for the test
		test_dis_gene_matrix_resp = test_responders_T[sort_index[:12]]
		test_dis_gene_matrix_non = test_non_responders_T[sort_index[:12]]
		test_both_matrix = np.concatenate((test_dis_gene_matrix_resp, test_dis_gene_matrix_non), axis=1)

		training_samp_names = training_resp_samp_names + training_non_resp_samp_names
		test_samp_names = test_resp_samp_names + test_non_resp_samp_names
		'''


		##########
		#HEATMAP
		'''
		plt.figure(1)

		
		cdict = {'red' : ((0., 1., 1.), (0.45, 0., 0.), (1., 0., 0.)),
				'green': ((0., 0., 0.), (0.55, 0., 0.), (1., 1., 1.)), 
				'blue' : ((0., 0., 0.), (0.9, 0., 0.), (1., 0., 0.))
				}
		
		#cdict = {'red' : ((0., 0., 0.), (0.45, 0., 0.), (1., 0., 0.)),
		#		'green': ((0., 0., 0.), (0.5, 0.5, 0.5), (1., 1., 1.)), 
		#		'blue' : ((0., 0., 0.), (0.9, 0., 0.), (1., 0., 0.))
		#		}

		my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,7)

		#list_of_samps = ['resp', 'resp', 'resp', 'resp', 'resp', 'resp', 'resp', 'resp', 
		#				'non_resp',	'non_resp',	'non_resp',	'non_resp',	'non_resp',	'non_resp',	'non_resp',	'non_resp',	]
		#make heat map of significant genes/samples

		#plt.subplot(121)
		plt.pcolor(both_matrix, cmap=my_cmap, alpha=0.8)
		plt.xticks(np.arange(0,len(both_matrix[0]))+0.5,training_samp_names, rotation='vertical', fontsize=5)
		plt.yticks(np.arange(0,12)+0.5, top_gene_names, fontsize=8)
		plt.tick_params(axis='both', which='both', left='off', right='off', bottom='off', top='off')
		plt.title('Expression Pattern of Distinguishing\nGenes in Training Data', fontsize='small')

		plt.colorbar(shrink=0.5)

		#list_of_samps = ['resp', 'resp', 'resp', 'resp', 'resp', 'resp', 'resp',  
		#			'non_resp',	'non_resp',	'non_resp',	'non_resp',	'non_resp',	'non_resp', 'non_resp']
		#plot test expressions to the right of the training
		
		plt.subplot(122)
		plt.pcolor(test_both_matrix, cmap=my_cmap, alpha=0.8)
		plt.xticks(np.arange(0,len(both_matrix[0]))+0.5,test_samp_names, rotation='vertical', fontsize=5)
		#plt.yticks(np.arange(0,12)+0.5, top_gene_names, fontsize=8)
		plt.tick_params(axis='both', which='both', left='off', right='off', bottom='off', top='off', labelleft='off')
		plt.title('Expression Pattern of Distinguishing\nGenes in Test Data', fontsize='small')
		



		#plt.annotate('Responders', (0.24,0.014), xycoords='figure fraction', fontsize=7)
		#plt.annotate('Non-Responders', (0.43,0.014), xycoords='figure fraction', fontsize=7)
		#plt.annotate('Responders', (0.63,0.014), xycoords='figure fraction', fontsize=7)
		#plt.annotate('Non-Responders', (0.76,0.014), xycoords='figure fraction', fontsize=7)

		plt.tight_layout()



		plt.savefig('heatmap_distinguish_genes.pdf')
		'''

		

		#make graph of CS graph


		'''
		#calculate prediciton score of each sample in training data and test data and plot like they did
		#use the top 12 genes
		training_responders_PS = []
		for sample in range(len(training_responders)):
			responder_votes = 0
			non_responder_votes = 0
			for gene_index in sort_index[:12]:
				#V
				vote_weight = abs(training_responders[sample][gene_index] - ((resp_mean[gene_index] + non_resp_mean[gene_index])/2))
				if abs(training_responders[sample][gene_index] - resp_mean[gene_index]) < abs(training_responders[sample][gene_index] - non_resp_mean[gene_index]):
					responder_votes += vote_weight
				else:
					non_responder_votes += vote_weight
			#PS
			PS = ((responder_votes - non_responder_votes) / (responder_votes + non_responder_votes)) * 100
			training_responders_PS.append(PS)

		#print
		#print training_responders_PS
		#print
		correct = []
		for guess in training_responders_PS:
			if guess > 0:
				correct.append(1)
			else:
				correct.append(0)


		training_non_responders_PS = []
		for sample in range(len(training_non_responders)):
			responder_votes = 0
			non_responder_votes = 0
			for gene_index in sort_index[:12]:
				#V
				vote_weight = abs(training_non_responders[sample][gene_index] - ((resp_mean[gene_index] + non_resp_mean[gene_index])/2))
				if abs(training_non_responders[sample][gene_index] - resp_mean[gene_index]) < abs(training_non_responders[sample][gene_index] - non_resp_mean[gene_index]):
					responder_votes += vote_weight
				else:
					non_responder_votes += vote_weight
			#PS
			PS = ((responder_votes - non_responder_votes) / (responder_votes + non_responder_votes)) * 100
			training_non_responders_PS.append(PS)

		#print training_non_responders_PS
		for guess in training_non_responders_PS:
			if guess < 0:
				correct.append(1)
			else:
				correct.append(0)


		print 'training accuracy = ' + str(np.mean(correct))



		test_responders_PS = []
		for sample in range(len(test_responders)):
			responder_votes = 0
			non_responder_votes = 0
			for gene_index in sort_index[:12]:
				#V
				vote_weight = abs(test_responders[sample][gene_index] - ((resp_mean[gene_index] + non_resp_mean[gene_index])/2))
				if abs(test_responders[sample][gene_index] - resp_mean[gene_index]) < abs(test_responders[sample][gene_index] - non_resp_mean[gene_index]):
					responder_votes += vote_weight
				else:
					non_responder_votes += vote_weight
			#PS
			PS = ((responder_votes - non_responder_votes) / (responder_votes + non_responder_votes)) * 100
			test_responders_PS.append(PS)

		correct = []
		for guess in test_responders_PS:
			if guess > 0:
				correct.append(1)
			else:
				correct.append(0)


		test_non_responders_PS = []
		for sample in range(len(test_non_responders)):
			responder_votes = 0
			non_responder_votes = 0
			for gene_index in sort_index[:12]:
				#V
				vote_weight = abs(test_non_responders[sample][gene_index] - ((resp_mean[gene_index] + non_resp_mean[gene_index])/2))
				if abs(test_non_responders[sample][gene_index] - resp_mean[gene_index]) < abs(test_non_responders[sample][gene_index] - non_resp_mean[gene_index]):
					responder_votes += vote_weight
				else:
					non_responder_votes += vote_weight
			#PS
			PS = ((responder_votes - non_responder_votes) / (responder_votes + non_responder_votes)) * 100
			test_non_responders_PS.append(PS)

		for guess in test_non_responders_PS:
			if guess < 0:
				correct.append(1)
			else:
				correct.append(0)


		print 'test accuracy = ' + str(np.mean(correct))
		'''

	def takata_fit(x, y, gene_names):

		#now genes will be rows and samples are columns
		x = x.T

		gene_set = ['SPRY1', 'OSBPL11', 'ZNF107', 'LIN7C', 'WDR90', 'SLC22A18', 'PNPO', 'CRKL']

		index_of_genes = []
		for gene in gene_set:
			index_of_genes.append(gene_names.index(gene))

		#get expressions for just those genes in each of the groups
		x = x[index_of_genes]

		#now genes will be columns and samples are rows
		x = x.T

		#calculate mean and std of training group for each gene
		resp_gene_expr = [[] for i in xrange(len(x[0]))]
		non_resp_gene_expr = [[] for i in xrange(len(x[0]))]

		for samp in range(len(x)):
			for gene in range(len(x[0])):
				if y[samp] == 0.0:
					resp_gene_expr[gene].append(x[samp][gene])
				else:
					non_resp_gene_expr[gene].append(x[samp][gene])

		mean_resp_gene_expr = []
		mean_non_resp_gene_expr = []

		for gene in range(len(x[0])):
			mean_resp_gene_expr.append(np.mean(resp_gene_expr[gene]))
			mean_non_resp_gene_expr.append(np.mean(non_resp_gene_expr[gene]))


		return mean_resp_gene_expr, mean_non_resp_gene_expr

	def takata_score(x, y, gene_names, mean_resp_gene_expr, mean_non_resp_gene_expr):
		
		#now genes will be rows and samples are columns
		x = x.T

		gene_set = ['SPRY1', 'OSBPL11', 'ZNF107', 'LIN7C', 'WDR90', 'SLC22A18', 'PNPO', 'CRKL']

		index_of_genes = []
		for gene in gene_set:
			index_of_genes.append(gene_names.index(gene))

		#get expressions for just those genes in each of the groups
		x = x[index_of_genes]

		#now genes will be columns and samples are rows
		x = x.T

		#calculate PS 
		PS_scores = []
		for sample in range(len(x)):
			responder_votes = 0
			non_responder_votes = 0
			for gene in range(len(x[0])):
				#V
				vote_weight = abs(x[sample][gene] - ((mean_resp_gene_expr[gene] + mean_non_resp_gene_expr[gene])/2))
				if abs(mean_resp_gene_expr[gene] - x[sample][gene] ) < abs(mean_non_resp_gene_expr[gene] - x[sample][gene]):
					responder_votes += vote_weight
				else:
					non_responder_votes += vote_weight
			#PS
			PS = ((responder_votes - non_responder_votes) / (responder_votes + non_responder_votes))
			PS_scores.append(PS)

		correct = []
		for i in range(len(PS_scores)):
			if (PS_scores[i] > 0 and y[i] == 0.0) or (PS_scores[i] < 0 and y[i] == 1.0):
				correct.append(1)
			else:
				correct.append(0)

		return np.mean(correct)

	def remove_low_expr_genes(training_responders, training_non_responders, test_responders, test_non_responders, gene_names):


		#cutoff = -0.05
		cutoff = -0.1
		genes_to_keep = []
		#make genes the rows
		training_responders_T = training_responders.T
		training_non_responders_T = training_non_responders.T
		test_responders_T = test_responders.T
		test_non_responders_T = test_non_responders.T
		print training_responders_T.shape
		print training_non_responders_T.shape
		#print training_data_T.shape
		for gene in range(len(training_responders_T)):
			above_cutoff = []
			for samp in range(len(training_responders_T[0])):
				if training_responders_T[gene][samp] > cutoff:
					above_cutoff.append(1)
				else:
					above_cutoff.append(0)
			if np.mean(above_cutoff) > 0.6 and np.std(training_responders_T[gene]) > 0.00001:
				genes_to_keep.append(gene)

		training_responders_T = training_responders_T[genes_to_keep]
		training_non_responders_T = training_non_responders_T[genes_to_keep]
		test_responders_T = test_responders_T[genes_to_keep]
		test_non_responders_T = test_non_responders_T[genes_to_keep]

		new_gene_names = []
		for i in genes_to_keep:
			new_gene_names.append(gene_names[i])
		gene_names = new_gene_names

		print training_responders_T.shape
		print training_non_responders_T.shape

		genes_to_keep = []
		for gene in range(len(training_non_responders_T)):
			above_cutoff = []
			for samp in range(len(training_non_responders_T[0])):
				if training_non_responders_T[gene][samp] > cutoff:
					above_cutoff.append(1)
				else:
					above_cutoff.append(0)
			if np.mean(above_cutoff) > 0.6 and np.std(training_non_responders_T[gene]) > 0.00001:
				genes_to_keep.append(gene)

		training_responders_T = training_responders_T[genes_to_keep]
		training_non_responders_T = training_non_responders_T[genes_to_keep]
		test_responders_T = test_responders_T[genes_to_keep]
		test_non_responders_T = test_non_responders_T[genes_to_keep]

		new_gene_names = []
		new_gene_names = []
		for i in genes_to_keep:
			new_gene_names.append(gene_names[i])
		gene_names = new_gene_names

		training_responders = training_responders_T.T
		training_non_responders = training_non_responders_T.T
		test_responders = test_responders_T.T
		test_non_responders = test_non_responders_T.T
		
		print training_responders_T.shape
		print training_non_responders_T.shape



		#this is where I want to remove genes with lwo standard deviation
		#training_data = np.concatenate((training_responders, training_non_responders), axis=0)	


		return training_responders, training_non_responders, test_responders, test_non_responders, gene_names

	def learn_and_test2(training_responders, training_non_responders, test_responders, test_non_responders, gene_names):


		training_data = np.concatenate((training_responders, training_non_responders), axis=0)
		targets_resp = [0]*len(training_responders)
		targets_non_resp = [1]*len(training_non_responders)
		training_targets = np.array(targets_resp+targets_non_resp).T

		test_data = np.concatenate((test_responders, test_non_responders), axis=0)
		targets_resp = [0]*len(test_responders)
		targets_non_resp = [1]*len(test_non_responders)
		test_targets = np.array(targets_resp+targets_non_resp).T

		
		#Logistic Regression
		clf = linear_model.LogisticRegression(penalty='l1')
		clf.fit(training_data, training_targets)
		pprint.pprint(clf.coef_)
		print 'Log Reg Training Score: ' + str(clf.score(training_data, training_targets))
		print 'Log Reg Test Score: ' + str(clf.score(test_data, test_targets))
		a = clf.predict(training_responders)
		b = clf.predict(training_non_responders)
		c = clf.predict(test_responders)
		d = clf.predict(test_non_responders)
		log_reg_training_predictions = list(a) + list(b)
		log_reg_test_predictions = list(c) + list(d)

		#Nearest Neighbors
		clf = KNeighborsClassifier(n_neighbors=3)
		clf.fit(training_data, training_targets)
		a = clf.predict(training_responders)
		b = clf.predict(training_non_responders)
		c = clf.predict(test_responders)
		d = clf.predict(test_non_responders)
		print 'Nearest Training Score: ' + str(clf.score(training_data, training_targets))
		print 'Nearest Test Score: ' + str(clf.score(test_data, test_targets))
		knn_training_predictions = list(a) + list(b)
		knn_test_predictions = list(c) + list(d)



		#Naive Bayes
		clf = GaussianNB()
		clf.fit(training_data, training_targets) 
		a = clf.predict(training_responders)
		b = clf.predict(training_non_responders)
		c = clf.predict(test_responders)
		d = clf.predict(test_non_responders)
		#print clf.predict_log_proba(training_data)
		print 'GNB Training Score: ' + str(clf.score(training_data, training_targets))
		#print clf.predict_log_proba(training_data)
		print 'GNB Training Score: ' + str(clf.score(test_data, test_targets))
		gnb_training_predictions = list(a) + list(b)
		gnb_test_predictions = list(c) + list(d)


		predictions = combine_predictions(log_reg_training_predictions, knn_training_predictions, gnb_training_predictions)
		print 'Combined Trainig score ' + str(score(predictions, training_targets))
		predictions = combine_predictions(log_reg_test_predictions, knn_test_predictions, gnb_test_predictions)
		print predictions
		print test_targets
		print 'Combined Test score ' + str(score(predictions, test_targets))

		'''
		#Linear Regression
		clf = linear_model.LinearRegression()
		clf.fit(training_data, training_targets)
		a = clf.predict(training_responders)
		b = clf.predict(training_non_responders)
		c = clf.predict(test_responders)
		d = clf.predict(test_non_responders)
		'''
		

		print
		clf = RandomForestClassifier()
		clf.fit(training_data, training_targets) 
		a = clf.predict(training_responders)
		b = clf.predict(training_non_responders)
		c = clf.predict(test_responders)
		d = clf.predict(test_non_responders)
		#print clf.predict_log_proba(training_data)
		print 'Random forest Training Score: ' + str(clf.score(training_data, training_targets))
		#print clf.predict_log_proba(training_data)
		print 'Random forest Training Score: ' + str(clf.score(test_data, test_targets))
		#gnb_training_predictions = list(a) + list(b)
		#gnb_test_predictions = list(c) + list(d)


		#clf = SVC()
		#clf.fit(training_data, training_targets)
		#print 'SVC Training Score: ' + str(clf.score(training_data, training_targets))
		#print 'SVC Test Score: ' + str(clf.score(test_data, test_targets))
		

		'''
		plt.figure(3)
		plt.scatter(np.ones(len(a)), a, c='blue')
		plt.scatter(np.ones(len(c))*2, c, c='lightblue')
		plt.scatter(np.ones(len(b))*4, b, c='red')	
		plt.scatter(np.ones(len(d))*3, d, c='pink')

		plt.title('Prediction Scores')
		plt.ylabel('Predictions')
		plt.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off')

		plt.annotate('Responders', (0.26,0.07), xycoords='figure fraction')
		plt.annotate('Non-Responders', (0.64,0.07), xycoords='figure fraction')
		plt.annotate('Training', (0.2,0.03), xycoords='figure fraction')
		plt.annotate('Test', (0.4,0.03), xycoords='figure fraction')
		plt.annotate('Test', (0.61,0.03), xycoords='figure fraction')
		plt.annotate('Training', (0.8,0.03), xycoords='figure fraction')

		plt.savefig('scatter_predictions_naivebayes.pdf')
		'''

	def pca_analysis(training_responders, training_non_responders):

		print 'PCA'

		print training_responders.shape
		print training_non_responders.shape

		training_data = np.concatenate((training_responders, training_non_responders), axis=0)

		pca = PCA(n_components=10)
		training_data = pca.fit_transform(training_data)

		training_responders = training_data[:len(training_responders)]
		training_non_responders = training_data[len(training_responders):]


		#print pca.components_
		print pca.explained_variance_ratio_
		print 'total variance explained ' + str(np.sum(pca.explained_variance_ratio_))

		print training_responders.shape
		print training_non_responders.shape


		#PLOT PCA
		#plt.figure(12)
		#plt.scatter(training_responders.T[0], training_responders.T[1], c='blue')
		#plt.scatter(training_non_responders.T[0], training_non_responders.T[1], c='red')
		#plt.savefig('scatter_pca.pdf')


		return training_responders, training_non_responders

	def most_predictive_genes(training_responders, training_non_responders, test_responders, test_non_responders, gene_names):

		print 'TOP GENES'
		#get median of each gene
		training_data = np.concatenate((training_responders, training_non_responders), axis=0)

		#median of the columns, which are the genes
		medians = np.median(training_data, axis=0)

		#find genes that has median that splits the classes best
		#since median splits samples in 2, if it splits one class well, it must split the other class well too
		score_positive_gene = []
		score_negative_gene = []
		for gene in range(len(training_responders.T)):
			score_positive_gene.append((training_responders.T[gene] > medians[gene]).sum())
			score_negative_gene.append((training_responders.T[gene] < medians[gene]).sum())

		len0 = len(training_responders)

		#score_positive_gene = np.array(score_positive_gene) / len0
		#score_negative_gene = np.array(score_negative_gene) / len0


		sort_index = np.argsort(score_positive_gene)
		sort_index2 = np.argsort(score_negative_gene)
		top_gene_names = []
		for i in range(12):
			top_gene_names.append(gene_names[sort_index[i]])
		print top_gene_names

		top_gene_names = []
		for i in range(12):
			top_gene_names.append(gene_names[sort_index2[i]])
		print top_gene_names

		a = sorted(score_positive_gene, reverse=True)
		b = sorted(score_negative_gene, reverse=True)
		print a[:20]
		print b[:20]

		#TODO, return these top genes data, so that I can train on them and get results. 


		genes_to_keep = []
		for i in range(len(score_positive_gene)):
			if score_positive_gene[i] >= ((len(training_data)/2) -1):
				genes_to_keep.append(i)
			elif score_negative_gene[i] >= ((len(training_data)/2) -1):
				genes_to_keep.append(i)

		training_responders = training_responders[:,genes_to_keep]
		training_non_responders = training_non_responders[:,genes_to_keep]
		test_responders = test_responders[:,genes_to_keep]
		test_non_responders = test_non_responders[:,genes_to_keep]
		gene_names2 = []
		for i in genes_to_keep:
			gene_names2.append(gene_names[i])
		gene_names = gene_names2


		#print training_responders.shape
		#print training_non_responders.shape
		#print gene_names

		return training_responders, training_non_responders, test_responders, test_non_responders, gene_names

	def combine_predictions(*argv):

		predictions = np.array(argv)
		
		combined_prediction = stats.mode(predictions)[0][0]

		return combined_prediction

	def score(prediction, training_targets):

		#print prediction

		count = 0
		#count2 = 0
		for i in range(len(prediction)):
			if prediction[i] == training_targets[i]:
				count += 1
			#if prediction[i] == test_targets[i]:
			#	count2 += 1


		training_score = float(count) / len(prediction) 
		#test_score = float(count2) / len(prediction)


		return training_score
		#print 'Score ' + str(training_score)
		#print 'Test score ' + str(test_score)

	

	def train_models_and_return_predictions(training_data, responses3, validation_sample, validation_response, gene_names):

		
		list_of_model_training_accuraccies = []
		list_of_model_validation_accuraccies = []

		predictions = []

		'''
		#Logistic Regression L1
		clf = linear_model.LogisticRegression(penalty='l1')
		clf.fit(training_data, responses3)
		training_score = clf.score(training_data, responses3)
		validation_score = clf.score([validation_sample], [validation_response])
		list_of_model_training_accuraccies.append(training_score)
		list_of_model_validation_accuraccies.append(validation_score)

		#Logistic Regression L2 C = 1
		clf = linear_model.LogisticRegression(penalty='l2', C=1.0)
		clf.fit(training_data, responses3)
		training_score = clf.score(training_data, responses3)
		validation_score = clf.score([validation_sample], [validation_response])
		list_of_model_training_accuraccies.append(training_score)
		list_of_model_validation_accuraccies.append(validation_score)

		#Logistic Regression L2 C = 0.5
		clf = linear_model.LogisticRegression(penalty='l2', C=0.5)
		clf.fit(training_data, responses3)
		training_score = clf.score(training_data, responses3)
		validation_score = clf.score([validation_sample], [validation_response])
		list_of_model_training_accuraccies.append(training_score)
		list_of_model_validation_accuraccies.append(validation_score)

		#Logistic Regression L2 C = 2
		clf = linear_model.LogisticRegression(penalty='l2', C=2.0)
		clf.fit(training_data, responses3)
		training_score = clf.score(training_data, responses3)
		validation_score = clf.score([validation_sample], [validation_response])
		list_of_model_training_accuraccies.append(training_score)
		list_of_model_validation_accuraccies.append(validation_score)

		'''

		#Naive Bayes
		clf = GaussianNB()
		clf.fit(training_data, responses3) 
		training_score = clf.score(training_data, responses3)
		validation_score = clf.score([validation_sample], [validation_response])
		list_of_model_training_accuraccies.append(training_score)
		list_of_model_validation_accuraccies.append(validation_score)
		prediction = clf.predict([validation_sample])
		predictions.append(prediction)

		'''
		#3 Nearest Neighbors
		clf = KNeighborsClassifier(n_neighbors=3)
		clf.fit(training_data, responses3)
		training_score = clf.score(training_data, responses3)
		validation_score = clf.score([validation_sample], [validation_response])
		list_of_model_training_accuraccies.append(training_score)
		list_of_model_validation_accuraccies.append(validation_score)

		#5 Nearest Neighbors
		clf = KNeighborsClassifier(n_neighbors=5)
		clf.fit(training_data, responses3)
		training_score = clf.score(training_data, responses3)
		validation_score = clf.score([validation_sample], [validation_response])
		list_of_model_training_accuraccies.append(training_score)
		list_of_model_validation_accuraccies.append(validation_score)

		#Random Forest
		clf = RandomForestClassifier()
		clf.fit(training_data, responses3) 
		training_score = clf.score(training_data, responses3)
		validation_score = clf.score([validation_sample], [validation_response])
		list_of_model_training_accuraccies.append(training_score)
		list_of_model_validation_accuraccies.append(validation_score)

		'''

		#SVC
		clf = SVC()
		clf.fit(training_data, responses3)
		training_score = clf.score(training_data, responses3)
		validation_score = clf.score([validation_sample], [validation_response])
		list_of_model_training_accuraccies.append(training_score)
		list_of_model_validation_accuraccies.append(validation_score)
		prediction = clf.predict([validation_sample])
		predictions.append(prediction)



		'''
		#GBC
		clf = GradientBoostingClassifier()
		clf.fit(training_data, responses3)
		training_score = clf.score(training_data, responses3)
		validation_score = clf.score([validation_sample], [validation_response])
		list_of_model_training_accuraccies.append(training_score)
		list_of_model_validation_accuraccies.append(validation_score)
		


		#Takata
		mean_resp_gene_expr, mean_non_resp_gene_expr = takata_fit(training_data, responses3, gene_names)
		training_score = takata_score(training_data, responses3, gene_names, mean_resp_gene_expr, mean_non_resp_gene_expr)
		validation_score = takata_score(np.array([validation_sample]), np.array([validation_response]), gene_names, mean_resp_gene_expr, mean_non_resp_gene_expr)
		list_of_model_training_accuraccies.append(training_score)
		list_of_model_validation_accuraccies.append(validation_score)
		'''

		predictions.append(validation_response)

		return list_of_model_training_accuraccies, list_of_model_validation_accuraccies, predictions

	def add_responses_to_data(data, responses):

		print 'data shape ' + str(data.shape)

		print 'responses shape ' + str(np.array([responses]).shape)

		data = np.concatenate((data, np.array([responses]).T), axis=1)

		print 'data shape ' + str(data.shape)


		return data

	def random_gene_set(data, gene_names, size):

		rand_set = random.sample(range(len(data[0])), size)
		rand_set.sort()
		#print rand_set

		data = data.T
		data = data[rand_set]
		data = data.T

		new_gene_names = []
		for i in rand_set:
			new_gene_names.append(gene_names[i])

		print 'Random gene set, data shape ' + str(data.shape)

		return data, new_gene_names

	def train_models(training_data, responses3, validation_sample, validation_response, gene_names):

		
		list_of_model_training_accuraccies = []
		list_of_model_validation_accuraccies = []
		list_of_times = []

		
		#Logistic Regression L1
		start = time()
		clf = linear_model.LogisticRegression(penalty='l1')
		clf.fit(training_data, responses3)
		training_score = clf.score(training_data, responses3)
		validation_score = clf.score([validation_sample], [validation_response])
		list_of_model_training_accuraccies.append(training_score)
		list_of_model_validation_accuraccies.append(validation_score)
		end = time()
		list_of_times.append(end - start)

		#Logistic Regression L2 C = 0.5
		start = time()
		clf = linear_model.LogisticRegression(penalty='l2', C=0.5)
		clf.fit(training_data, responses3)
		training_score = clf.score(training_data, responses3)
		validation_score = clf.score([validation_sample], [validation_response])
		list_of_model_training_accuraccies.append(training_score)
		list_of_model_validation_accuraccies.append(validation_score)
		end = time()
		list_of_times.append(end - start)


		#Logistic Regression L2 C = 1
		start = time()
		clf = linear_model.LogisticRegression(penalty='l2', C=1.0)
		clf.fit(training_data, responses3)
		training_score = clf.score(training_data, responses3)
		validation_score = clf.score([validation_sample], [validation_response])
		list_of_model_training_accuraccies.append(training_score)
		list_of_model_validation_accuraccies.append(validation_score)
		end = time()
		list_of_times.append(end - start)

		#Logistic Regression L2 C = 2
		start = time()
		clf = linear_model.LogisticRegression(penalty='l2', C=2.0)
		clf.fit(training_data, responses3)
		training_score = clf.score(training_data, responses3)
		validation_score = clf.score([validation_sample], [validation_response])
		list_of_model_training_accuraccies.append(training_score)
		list_of_model_validation_accuraccies.append(validation_score)
		end = time()
		list_of_times.append(end - start)

		#SGD classifier with elastic net regularization and hinge loss
		start = time()
		clf = linear_model.SGDClassifier(penalty='elasticnet')
		clf.fit(training_data, responses3)
		training_score = clf.score(training_data, responses3)
		validation_score = clf.score([validation_sample], [validation_response])
		list_of_model_training_accuraccies.append(training_score)
		list_of_model_validation_accuraccies.append(validation_score)
		end = time()
		list_of_times.append(end - start)


		#Naive Bayes
		start = time()
		clf = GaussianNB()
		clf.fit(training_data, responses3) 
		training_score = clf.score(training_data, responses3)
		validation_score = clf.score([validation_sample], [validation_response])
		list_of_model_training_accuraccies.append(training_score)
		list_of_model_validation_accuraccies.append(validation_score)
		end = time()
		list_of_times.append(end - start)

		
		#3 Nearest Neighbors
		start = time()
		clf = KNeighborsClassifier(n_neighbors=3)
		clf.fit(training_data, responses3)
		training_score = clf.score(training_data, responses3)
		validation_score = clf.score([validation_sample], [validation_response])
		list_of_model_training_accuraccies.append(training_score)
		list_of_model_validation_accuraccies.append(validation_score)
		end = time()
		list_of_times.append(end - start)

		#5 Nearest Neighbors
		start = time()
		clf = KNeighborsClassifier(n_neighbors=5)
		clf.fit(training_data, responses3)
		training_score = clf.score(training_data, responses3)
		validation_score = clf.score([validation_sample], [validation_response])
		list_of_model_training_accuraccies.append(training_score)
		list_of_model_validation_accuraccies.append(validation_score)
		end = time()
		list_of_times.append(end - start)

		#Random Forest
		start = time()
		clf = RandomForestClassifier()
		clf.fit(training_data, responses3) 
		training_score = clf.score(training_data, responses3)
		validation_score = clf.score([validation_sample], [validation_response])
		list_of_model_training_accuraccies.append(training_score)
		list_of_model_validation_accuraccies.append(validation_score)
		end = time()
		list_of_times.append(end - start)

		#SVC
		start = time()
		clf = SVC()
		clf.fit(training_data, responses3)
		training_score = clf.score(training_data, responses3)
		validation_score = clf.score([validation_sample], [validation_response])
		list_of_model_training_accuraccies.append(training_score)
		list_of_model_validation_accuraccies.append(validation_score)


		#LDA
		start = time()
		clf = lda.LDA()
		clf.fit(training_data, responses3)
		training_score = clf.score(training_data, responses3)
		validation_score = clf.score([validation_sample], [validation_response])
		list_of_model_training_accuraccies.append(training_score)
		list_of_model_validation_accuraccies.append(validation_score)
		end = time()
		list_of_times.append(end - start)

		'''
		#QDA
		start = time()
		clf = qda.QDA()
		clf.fit(training_data, responses3)
		training_score = clf.score(training_data, responses3)
		validation_score = clf.score([validation_sample], [validation_response])
		list_of_model_training_accuraccies.append(training_score)
		list_of_model_validation_accuraccies.append(validation_score)
		end = time()
		list_of_times.append(end - start)
		'''

		
		#GBC
		start = time()
		clf = GradientBoostingClassifier()
		clf.fit(training_data, responses3)
		training_score = clf.score(training_data, responses3)
		validation_score = clf.score([validation_sample], [validation_response])
		list_of_model_training_accuraccies.append(training_score)
		list_of_model_validation_accuraccies.append(validation_score)
		end = time()
		list_of_times.append(end - start)
		

		'''
		#Takata
		start = time()
		mean_resp_gene_expr, mean_non_resp_gene_expr = takata_fit(training_data, responses3, gene_names)
		training_score = takata_score(training_data, responses3, gene_names, mean_resp_gene_expr, mean_non_resp_gene_expr)
		validation_score = takata_score(np.array([validation_sample]), np.array([validation_response]), gene_names, mean_resp_gene_expr, mean_non_resp_gene_expr)
		list_of_model_training_accuraccies.append(training_score)
		list_of_model_validation_accuraccies.append(validation_score)
		end = time()
		list_of_times.append(end - start)
		'''


		return list_of_model_training_accuraccies, list_of_model_validation_accuraccies, list_of_times

	def remove_low_std_genes(cutoff, training_responders, training_non_responders, test_responders, test_non_responders, gene_names):

		print 'Training responders shape before std filter of : ' + str(cutoff) + ' ' + str(training_responders.shape)

		training_data = np.concatenate((training_responders,training_non_responders), axis=0)
		training_data_T = training_data.T

		genes_to_keep = []
		new_gene_names = []
		for gene in range(len(training_data_T)):
			if np.std(training_data_T[gene]) > cutoff:
				genes_to_keep.append(gene)
				new_gene_names.append(gene_names[gene])

		training_data_T = training_data_T[genes_to_keep]
		training_data = training_data_T.T
		training_responders = training_data[:len(training_responders)]
		training_non_responders = training_data[len(training_responders):]
		test_responders = test_responders.T[genes_to_keep].T
		test_non_responders = test_non_responders.T[genes_to_keep].T

		print 'Training responders shape after std filter: ' + str(training_responders.shape)


		return training_responders, training_non_responders, test_responders, test_non_responders, new_gene_names

	def remove_low_mean_genes(cutoff, training_responders, training_non_responders, test_responders, test_non_responders, gene_names):

		print 'Training responders shape before mean filter of : ' + str(cutoff) + ' ' + str(training_responders.shape)

		training_data = np.concatenate((training_responders,training_non_responders), axis=0)
		training_data_T = training_data.T

		genes_to_keep = []
		new_gene_names = []
		for gene in range(len(training_data_T)):
			if np.mean(training_data_T[gene]) > cutoff:
				genes_to_keep.append(gene)
				new_gene_names.append(gene_names[gene])

		training_data_T = training_data_T[genes_to_keep]
		training_data = training_data_T.T
		training_responders = training_data[:len(training_responders)]
		training_non_responders = training_data[len(training_responders):]
		test_responders = test_responders.T[genes_to_keep].T
		test_non_responders = test_non_responders.T[genes_to_keep].T

		print 'Training responders shape after mean filter: ' + str(training_responders.shape)


		return training_responders, training_non_responders, test_responders, test_non_responders, new_gene_names



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

	reponse = np.array(response)

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

	#print 'Numb of genes ' + str(len(data))
	#print 'Numb of samples ' + str(len(data[0]))


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


	return data, gene_names, valid_samp


def neural_net(data, responses, early_stopping_samples, early_stopping_responses, validation_sample, validation_response):

		
		#input is a list of tuples (x, y), x=array of input data, y=expected output
		training_data = []
		for i in range(len(data)):
			training_data.append((np.reshape(data[i], (-1,1)), np.reshape(np.array([responses[i]]), (-1,1))))


		evaluation_data = []
		for i in range(len(early_stopping_samples)):
			evaluation_data.append((np.reshape(early_stopping_samples[i], (-1,1)), np.reshape(np.array([early_stopping_responses[i]]), (-1,1))))

		
		epochs=500
		mini_batch_size=30
		learn_rate=0.05
		lmbda=0.01
		first_layer=10
		second_layer=10
		third_layer=100
		fourth_layer=3
		fifth_layer=3

		numb_of_non_resp = 0.0
		for i in responses:
			if i == 1.0:
				numb_of_non_resp +=1.0
		decision_boundary = numb_of_non_resp / float(len(responses))
		decision_boundary = 0.5

		#hyperparameters = [epochs,mini_batch_size,learn_rate,lmbda,first_layer,second_layer]


		net = NN_cc.Network([len(training_data[0][0]), first_layer, second_layer, 1])
		evaluation_cost, evaluation_accuracy, training_cost, training_accuracy, j = net.SGD(training_data,
																						evaluation_data=evaluation_data,
																						epochs=epochs, 
																						mini_batch_size=mini_batch_size, 
																						learn_rate =learn_rate, 
																						lmbda = lmbda,
																						monitor_evaluation_cost=True,
																						monitor_training_cost=True,
																						boundary=decision_boundary)



		score = net.accuracy([(np.reshape(validation_sample, (-1,1)), np.reshape(np.array(validation_response), (-1,1)))], decision_boundary)

		#if score == 1.0:
		#	print 'Got it right'
		#else:
		#	print 'Got it wrong'

		#print

		return score, j, training_cost, evaluation_cost


def doWork(valid_sample):

	#f = open('progress/' + str(valid_sample),'w')
	#f.write("Validation sample : " + str(valid_sample) + ' beginning...') 
	#f.close()
	print "Validation sample : " + str(valid_sample) + ' beginning...'

	#print 'Assembling data...'

	ids, responses, samps, first_row = assemble_data('rawdata_main_genes_only.csv', 'ClinpathDataBernCohort.xlsx', 'SampleTranslation.xlsx')
	data, sample_names = get_all_expressions('rawdata_main_genes_only.csv', first_row, samps, ids, responses)
	gene_names = get_gene_names('rawdata_main_genes_only.csv')


	numb_of_samples = 52
	numb_of_features = 9000
	numb_of_sig_features = 0

	#data, responses = make_simulated(numb_of_samples, numb_of_features, numb_of_sig_features, 15)



	#print 'Data assembly complete.'

	#make a new list of the reamaining samples -> training samples
	training_data = np.delete(data, valid_sample, 0)
	responses3 = np.delete(responses, valid_sample, 0)

	training_data, new_gene_names, valid_samp = filter_genes2(0.1, 0.1, training_data, gene_names, data[valid_sample])

	#get some early stopping samples
	#this should be random
	#TODO
	responder_early_stoppeers = []
	non_responder_early_stoppeers = []
	leftover = []
	for i in range(len(responses3)):
		if responses3[i] == 0.0 and len(responder_early_stoppeers) < 5:
			responder_early_stoppeers.append(i)
		elif responses3[i] == 1.0 and len(non_responder_early_stoppeers) < 5:
			non_responder_early_stoppeers.append(i)
		else:
			leftover.append(i)

	early_stopping_indices = responder_early_stoppeers + non_responder_early_stoppeers
	early_stopping_samples = training_data[early_stopping_indices]
	early_stopping_responses = responses3[early_stopping_indices]

	training_data = training_data[leftover]
	responses3 = responses3[leftover]

	#score is if it got validation sample correct
	score, epochs, training_cost, evaluation_cost = neural_net(training_data, responses3, early_stopping_samples, early_stopping_responses, valid_samp, responses[valid_sample])

	#for one of the validation runs, make a graph of training and validaiton error.
	if valid_sample == 1:

		plt.figure(56)
		x = range(len(training_cost))
		plt.plot(x, training_cost, ls='dashed')
		plt.plot(x, evaluation_cost)
		plt.xlabel('Epoch')
		plt.ylabel('Error')
		plt.savefig('NN_costs_simulated.pdf')



	print 'sample : ' + str(valid_sample) + ' complete.'

	#f = open('progress/' + str(valid_sample),'a')
	#f.write('sample : ' + str(valid_sample) + ' complete.')
	#f.close() 

	return score






if __name__ == "__main__":


	start = time()

	pool = Pool(processes=52)
	validation_samples = range(52)
	results = pool.map(doWork, validation_samples)

	end = time()
	total = (end - start)

	print '\nTotal time = ' + str(total)
	print 'Mean accuracy' + str(np.mean(results))


	f = open('results','a')
	f.write('\n----------------------\n')
	f.write('First layer = 20\n')
	f.write('Second layer = 20\n')
	#f.write('Third layer = 100\n')
	f.write('Learn rate = 0.05\n')
	f.write('Lmbda = 0.1\n')
	f.write('Max Epochs = 100\n')
	f.write('Mini batch size = 30\n')
	#f.write('Early stop = 100 inc epochs\n')
	#f.write('With new boundary\n')
	f.write('Run time = ' + str(total) + '\n')
	f.write('Accuracy = ' + str(np.mean(results)) + '\n')
	f.write('----------------------\n')
	f.close()