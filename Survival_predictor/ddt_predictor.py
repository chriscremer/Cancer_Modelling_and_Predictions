

import csv
import pprint
import random
from sklearn import linear_model
from sklearn import svm
from sklearn import neighbors
from sklearn import tree
from sklearn import ensemble
from sklearn import preprocessing

import NN_cc_regression_tanh

import numpy as np

clinical_file = '/data1/morrislab/ccremer/TCGA_data/rnaseqv2_clinical_onlytumours_25nov2014/Clinical/Biotab/nationwidechildrens.org_clinical_' 

f29 = open('update.txt', 'w')
f29.close()



###########################################################################
#
#  Gather Data
#
###########################################################################

print 'Gathering data...'


drug_list = [['cisplatin'], 
	['gemcitabine'],
	['bcg'],
	['paclitaxel', 'taxol', 'pacitaxol'],
	['methotrexate'],
	['platinum'],
	['vinorelbine'], 
	['halichondrin b', 'halichondrin'], 
	['mesna'], 
	['gemzar'], 
	['carboplatin'], 
	['docetaxel', 'docetaxel +/- zactima'], 
	['etoposide'], 
	['temsirolimus'],
	['vinblastine'],
	['adriamycin (doxorubicin)', 'doxorubicin', 'adriamycin'],
	['ifosfamide']]

#0 - ddt : death days to , This is actually the output that Im trying to predict
#1 - birth days to initial diagnosis
#2 -gender :#0 = female , 1 = male

#3 - cisplatin
#4 - gemcitabine
#5 - bcg
#6 - paclitaxel
#7 - methotrexate
#8 - platinum
#9 - vinorelbine
#10 - halichondrin b
#11 - mesna
#12 - gemzar
#13 - carboplatin
#14 - docetaxel
#15 - etoposide
#16 - temsirolimus
#17 - vinblastine
#18 - adriamycin (doxorubicin)
#19 - ifosfamide


variables = [[] for i in range(3)] #from patient file
drug_features = [[] for i in range(len(drug_list))]

#the columns where those features are found in the patient file
columns = [0]*3

patient_list = []

with open(clinical_file + 'patient_blca.txt', 'rU') as f:
	reader = csv.reader(f, delimiter='\t')
	count = 0
	for row in reader:

		#print count
		#count += 1

		#skip first two rows
		if row[0] == 'bcr_patient_barcode':
			#find index/column of features I want
			for column in range(len(row)):
				if row[column] == 'death_days_to':
					columns[0] = column
				#elif row[column] == 'age_at_diagnosis':
				elif row[column] == 'birth_days_to':				
					columns[1] = column
				elif row[column] == 'gender':
					columns[2] = column
				#elif row[column] == 'height_cm_at_diagnosis':
				#	columns[3] = column
			continue

		elif 'CDE_ID' in row[0]:
			continue
		

		#see if this patient/row has all the necessary data
		good_row = 1
		for column in columns:
			if row[column] == '[Not Available]' or row[column] == '[Not Applicable]' or row[column] == '[Unknown]':
				good_row = 0
				continue
		

		#if this patient has the necessary data, save the data
		if good_row == 1:
			patient_list.append(row[0])

			variables[0].append(float(row[columns[0]])) #days till death
			variables[1].append(float(row[columns[1]])*-1) #age , its in negative days
			if row[columns[2]] == 'MALE':
				variables[2].append(1.0)
			else:
				variables[2].append(0.0)
			#variables[3].append(float(row[columns[3]])) #height
		else:
			continue

		########################
		#get the drugs that the patient took
		########################

		patient_drugs = []
		patient_barcode = row[0]
		drug_name_column = 5
		with open(clinical_file + 'drug_blca.txt', 'rU') as f2:
			reader2 = csv.reader(f2, delimiter='\t')
			row_number = 0
			for row2 in reader2:
				if row2[0] == patient_barcode:
					if not row2[drug_name_column].lower() == '[not available]':
						patient_drugs.append(row2[drug_name_column].lower())
		#set all of the drug features to 0 for this patient, then change to 1 if they have it
		for list_of_drug in drug_features:
			list_of_drug.append(0)
		for drug_name in patient_drugs:
			drug_index = 0
			assigned = 0
			for specific_drug_name_list in drug_list:
				if drug_name in specific_drug_name_list:
					drug_features[drug_index][len(drug_features[drug_index])-1] = 1 
					assigned = 1
					continue
				drug_index += 1
			if assigned == 0:
				print drug_name

########################
#Get RNA expression data
########################

#make a dict of the mapping between patient barcode and rnaseq file name
patient_rnaseq_file_dict = {}
with open('/data1/morrislab/ccremer/TCGA_data/rnaseqv2_clinical_onlytumours_25nov2014/FILE_SAMPLE_MAP.txt', 'rU') as f:
	reader = csv.reader(f, delimiter='\t')
	count = 0
	for row in reader:
		if len(row) > 0:
			if '.genes.normalized_results' in row[0]:
				patient_rnaseq_file_dict[row[1]] = row[0]

#make a list of the rna seq files that match the patients that I'm using
patient_rnaseq_files = []
for patient_barcode in patient_list:
	assigned = 0
	for long_patient_barcode in patient_rnaseq_file_dict:
		if patient_barcode in long_patient_barcode:
			patient_rnaseq_files.append(patient_rnaseq_file_dict[long_patient_barcode])
			assigned = 1
			break
	if assigned == 0:
		print 'PROBLEM HERE'
		print patient_barcode



#list of all the gene ids
gene_ids = []
file1 = patient_rnaseq_files[0]
with open('/data1/morrislab/ccremer/TCGA_data/rnaseqv2_clinical_onlytumours_25nov2014/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/' + file1) as f:
	row_numb = 0
	for row in f:
		if row_numb == 0:
			row_numb+=1
			continue
		tmp = row.split()
		gene_ids.append(tmp[0])


rna_features = [[] for i in range(len(gene_ids))]


for file1 in patient_rnaseq_files:
	with open('/data1/morrislab/ccremer/TCGA_data/rnaseqv2_clinical_onlytumours_25nov2014/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/' + file1) as f:
		row_numb = 0
		gene_numb = 0
		for row in f:
			if row_numb == 0:
				row_numb+=1
				continue
			tmp = row.split()
			rna_features[gene_numb].append(float(tmp[1]))
			gene_numb += 1





#print len(rna_features)
#print len(rna_features[0])

#dimentionality reduction
#remove genes with low expression (max and avg) and variation

#avg_gene_expr = []
#std_gene_expr = []
#for gene_expr_list in rna_features:
#	avg_gene_expr.append(np.mean(gene_expr_list))
#	std_gene_expr.append(np.std(gene_expr_list))


my_randoms = random.sample(xrange(100), len(variables[0]))
#print my_randoms

#only height
#variables.pop(1)
#variables.pop(2) 

#remove age
#variables.pop(1)

#remove gender
#ariables.pop(2)

#assemble features
# so its features are the rows and sample/patient are the columns atm
features = variables[1:] + drug_features + rna_features #+ [my_randoms] # 
outputs = variables[0]


print 'numb 0f features ' + str(len(features))
print 'numb of samples ' + str(len(features[0]))
print 'Done.'


print 'Rescaling...'

features_matrix = np.array(features)
#transpose so samples are the rows
features_matrix = np.transpose(features_matrix)
#print 'shape'
#print features_matrix.shape

#scale to (0,1)
min_max_scaler = preprocessing.MinMaxScaler()
scaled_features = min_max_scaler.fit_transform(features_matrix)

#transpose so features are the rows again
scaled_features = np.transpose(scaled_features)
#print 'shape'
#print scaled_features.shape
outputs_array = np.array(outputs, ndmin=2)

print 'Done.'



print 'avg days ' + str(np.mean(outputs_array))
print 'std days ' + str(np.std(outputs_array))
print 'max days ' + str(np.amax(outputs_array))
print 'min days ' + str(np.amin(outputs_array))

print 'avg days ' + str(np.mean(scaled_features[0]))
print 'std days ' + str(np.std(scaled_features[0]))
print 'max days ' + str(np.amax(scaled_features[0]))
print 'min days ' + str(np.amin(scaled_features[0]))

#print outputs_array
#print scaled_features[0]
#print scaled_features[1]

#print outputs_array.shape
#print scaled_features.shape

output_and_features = np.concatenate((outputs_array, scaled_features))

#print output_and_features
print 'shape'
print output_and_features.shape


output_and_features_list = output_and_features.tolist()

avg_life = np.mean(outputs_array)



###########################################################################
#
#  Make different sets for cross-validation  --- should test on the train on every iter
#
###########################################################################

print 'Training and Testing Models...'

predictionIsAVG = []
least_squares_scores = []
svr_scores = []
neighbors_scores = []
tree_scores = []
forests_scores = []
neural_net_scores = []

predictionIsAVG_trainingSet = []
least_squares_scores_trainingSet = []
svr_scores_trainingSet = []
neighbors_scores_trainingSet = []
tree_scores_trainingSet = []
forests_scores_trainingSet = []
neural_net_scores_trainingSet = []


for fold_iteration in range(50):
	print 'Iteration: ' + str(fold_iteration)

	f29 = open('update.txt', 'a')
	f29.write('Iteration: ' + str(fold_iteration))	
	f29.close()
	
	#copying the features so that features doesnt change, set_b will
	#set_b = features[:]
	set_b = []
	for feat in output_and_features_list:
		set_b.append(list(feat))

	set_a = [[] for i in range(len(set_b))]

	#2 fold cross validation
	#so select half the samples to be used for vaidation
	for j in range(len(set_b[0])/2):
		rando = random.randint(0, len(set_b[0])-1)
		for feat in range(len(set_b)):
			set_a[feat].append(set_b[feat].pop(rando))

	test_x = np.array(set_a[1:]).transpose()
	test_y = np.array(set_a[0])

	training_x = np.array(set_b[1:]).transpose()
	training_y = np.array(set_b[0])

	#print 'Training data: ' + str(training_x.shape) + ' ' + str(training_y.shape)
	#print 'Test data: ' + str(test_x.shape) + ' ' + str(test_y.shape)


	###########################################################################
	#
	#  Fit and Test Models
	#
	###########################################################################


	for sample in range(len(test_x)):
		y = avg_life
		predictionIsAVG.append(abs(y- test_y[sample]))
	for sample in range(len(training_x)):
		y = avg_life
		predictionIsAVG_trainingSet.append(abs(y- training_y[sample]))


	least_squares = linear_model.LinearRegression()
	least_squares.fit(training_x, training_y)
	#print 'r squared ' + str(least_squares.score(test_x, test_y))
	for sample in range(len(test_x)):
		y = least_squares.predict(test_x[sample])
		least_squares_scores.append(abs(y- test_y[sample]))
	for sample in range(len(training_x)):
		y = least_squares.predict(training_x[sample])
		least_squares_scores_trainingSet.append(abs(y- training_y[sample]))


	svr = svm.SVR()
	svr.fit(training_x, training_y) 
	#print 'r squared ' + str(svr.score(test_x, test_y))
	for sample in range(len(test_x)):
		y = svr.predict(test_x[sample])
		svr_scores.append(abs(y- test_y[sample]))
	for sample in range(len(training_x)):
		y = svr.predict(training_x[sample])
		svr_scores_trainingSet.append(abs(y- training_y[sample]))


	neighrest_n = neighbors.KNeighborsRegressor()
	neighrest_n.fit(training_x, training_y) 
	#print 'r squared ' + str(neighrest_n.score(test_x, test_y))
	for sample in range(len(test_x)):
		y = neighrest_n.predict(test_x[sample])
		neighbors_scores.append(abs(y- test_y[sample]))
	for sample in range(len(training_x)):
		y = neighrest_n.predict(training_x[sample])
		neighbors_scores_trainingSet.append(abs(y- training_y[sample]))


	decision_tree = tree.DecisionTreeRegressor()
	decision_tree.fit(training_x, training_y) 
	#print 'r squared ' + str(decision_tree.score(test_x, test_y))
	for sample in range(len(test_x)):
		y = decision_tree.predict(test_x[sample])
		tree_scores.append(abs(y- test_y[sample]))
	for sample in range(len(training_x)):
		y = decision_tree.predict(training_x[sample])
		tree_scores_trainingSet.append(abs(y- training_y[sample]))


	random_forests = ensemble.RandomForestRegressor()
	random_forests.fit(training_x, training_y) 
	#print 'r squared ' + str(random_forests.score(test_x, test_y))
	for sample in range(len(test_x)):
		y = random_forests.predict(test_x[sample])
		forests_scores.append(abs(y- test_y[sample]))
	for sample in range(len(training_x)):
		y = random_forests.predict(training_x[sample])
		forests_scores_trainingSet.append(abs(y- training_y[sample]))

	###########################################################################
	# Switch the the two sets
	###########################################################################

	test_x = np.array(set_b[1:]).transpose()
	test_y = np.array(set_b[0])

	training_x = np.array(set_a[1:]).transpose()
	training_y = np.array(set_a[0])

	#print 'Training data: ' + str(training_x.shape) + ' ' + str(training_y.shape)
	#print 'Test data: ' + str(test_x.shape) + ' ' + str(test_y.shape)


	###########################################################################
	#
	#  Fit and Test Models
	#
	###########################################################################

	for sample in range(len(test_x)):
		y = avg_life
		predictionIsAVG.append(abs(y- test_y[sample]))
	for sample in range(len(training_x)):
		y = avg_life
		predictionIsAVG_trainingSet.append(abs(y- training_y[sample]))


	least_squares = linear_model.LinearRegression()
	least_squares.fit(training_x, training_y)
	#print 'r squared ' + str(least_squares.score(test_x, test_y))
	for sample in range(len(test_x)):
		y = least_squares.predict(test_x[sample])
		least_squares_scores.append(abs(y- test_y[sample]))
	for sample in range(len(training_x)):
		y = least_squares.predict(training_x[sample])
		least_squares_scores_trainingSet.append(abs(y- training_y[sample]))


	svr = svm.SVR()
	svr.fit(training_x, training_y) 
	#print 'r squared ' + str(svr.score(test_x, test_y))
	for sample in range(len(test_x)):
		y = svr.predict(test_x[sample])
		svr_scores.append(abs(y- test_y[sample]))
	for sample in range(len(training_x)):
		y = svr.predict(training_x[sample])
		svr_scores_trainingSet.append(abs(y- training_y[sample]))


	neighrest_n = neighbors.KNeighborsRegressor()
	neighrest_n.fit(training_x, training_y) 
	#print 'r squared ' + str(neighrest_n.score(test_x, test_y))
	for sample in range(len(test_x)):
		y = neighrest_n.predict(test_x[sample])
		neighbors_scores.append(abs(y- test_y[sample]))
	for sample in range(len(training_x)):
		y = neighrest_n.predict(training_x[sample])
		neighbors_scores_trainingSet.append(abs(y- training_y[sample]))


	decision_tree = tree.DecisionTreeRegressor()
	decision_tree.fit(training_x, training_y) 
	#print 'r squared ' + str(decision_tree.score(test_x, test_y))
	for sample in range(len(test_x)):
		y = decision_tree.predict(test_x[sample])
		tree_scores.append(abs(y- test_y[sample]))
	for sample in range(len(training_x)):
		y = decision_tree.predict(training_x[sample])
		tree_scores_trainingSet.append(abs(y- training_y[sample]))


	random_forests = ensemble.RandomForestRegressor()
	random_forests.fit(training_x, training_y) 
	#print 'r squared ' + str(random_forests.score(test_x, test_y))
	for sample in range(len(test_x)):
		y = random_forests.predict(test_x[sample])
		forests_scores.append(abs(y- test_y[sample]))
	for sample in range(len(training_x)):
		y = random_forests.predict(training_x[sample])
		forests_scores_trainingSet.append(abs(y- training_y[sample]))

	#training_data = [(x.reshape((len(training_x[0]), 1)), np.array([y]).reshape(1,1)) for x,y in zip(training_x, training_y)]
	#test_data = [(x.reshape((len(test_x[0]), 1)), np.array([y]).reshape(1,1)) for x,y in zip(test_x, test_y)]

	'''
	net = NN_cc_regression_tanh.Network([len(training_data[0][0]), 500, 500, 1])
	evaluation_cost, evaluation_accuracy, training_cost, training_accuracy = net.SGD(training_data, 10, 10, 0.001, lmbda=0.001, momemtum=0.01, evaluation_data=test_data, monitor_training_cost=True, monitor_evaluation_cost=True)
	for sample in range(len(test_x)):
		y = net.feedforward(test_x[sample])
		neural_net_scores.append(abs(y- test_y[sample]))
	'''



print '\nSupport Vector Regression'
#print svr_scores
print 'avg ' + str(np.mean(svr_scores))
print 'std ' + str(np.std(svr_scores))
print 'max ' + str(float(max(svr_scores)))
print 'min ' + str(float(min(svr_scores)))
print '\nSupport Vector Regression training'
#print svr_scores
print 'avg ' + str(np.mean(svr_scores_trainingSet))
print 'std ' + str(np.std(svr_scores_trainingSet))
print 'max ' + str(float(max(svr_scores_trainingSet)))
print 'min ' + str(float(min(svr_scores_trainingSet)))


print '\nNeighrest Neighbors'
print 'avg ' + str(np.mean(neighbors_scores))
print 'std ' + str(np.std(neighbors_scores))
print 'max ' + str(float(max(neighbors_scores)))
print 'min ' + str(float(min(neighbors_scores)))
print '\nNeighrest Neighbors training'
print 'avg ' + str(np.mean(neighbors_scores_trainingSet))
print 'std ' + str(np.std(neighbors_scores_trainingSet))
print 'max ' + str(float(max(neighbors_scores_trainingSet)))
print 'min ' + str(float(min(neighbors_scores_trainingSet)))



print '\nRandom Forests'
print 'avg ' + str(np.mean(forests_scores))
print 'std ' + str(np.std(forests_scores))
print 'max ' + str(float(max(forests_scores)))
print 'min ' + str(float(min(forests_scores)))
print '\nRandom Forests training'
print 'avg ' + str(np.mean(forests_scores_trainingSet))
print 'std ' + str(np.std(forests_scores_trainingSet))
print 'max ' + str(float(max(forests_scores_trainingSet)))
print 'min ' + str(float(min(forests_scores_trainingSet)))




print '\nDecision Tree'
print 'avg ' + str(np.mean(tree_scores))
print 'std ' + str(np.std(tree_scores))
print 'max ' + str(float(max(tree_scores)))
print 'min ' + str(float(min(tree_scores)))
print '\nDecision Tree training'
print 'avg ' + str(np.mean(tree_scores_trainingSet))
print 'std ' + str(np.std(tree_scores_trainingSet))
print 'max ' + str(float(max(tree_scores_trainingSet)))
print 'min ' + str(float(min(tree_scores_trainingSet)))



print '\nOrdinary Least Squares'
#print least_squares_scores
print 'avg ' + str(np.mean(least_squares_scores))
print 'std ' + str(np.std(least_squares_scores))
print 'max ' + str(max(least_squares_scores))
print 'min ' + str(min(least_squares_scores))
print '\nOrdinary Least Squares training'
#print least_squares_scores
print 'avg ' + str(np.mean(least_squares_scores_trainingSet))
print 'std ' + str(np.std(least_squares_scores_trainingSet))
print 'max ' + str(max(least_squares_scores_trainingSet))
print 'min ' + str(min(least_squares_scores_trainingSet))


print '\nLife Expectency = Average'
#print least_squares_scores
print 'avg ' + str(np.mean(predictionIsAVG))
print 'std ' + str(np.std(predictionIsAVG))
print 'max ' + str(max(predictionIsAVG))
print 'min ' + str(min(predictionIsAVG))
print '\nLife Expectency = Average training'
#print least_squares_scores
print 'avg ' + str(np.mean(predictionIsAVG_trainingSet))
print 'std ' + str(np.std(predictionIsAVG_trainingSet))
print 'max ' + str(max(predictionIsAVG_trainingSet))
print 'min ' + str(min(predictionIsAVG_trainingSet))


#print '\nNeural Net'
#print 'avg days off from real death ' + str(np.mean(neural_net_scores))

