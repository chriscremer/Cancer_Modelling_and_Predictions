


import numpy as np



def split_train_test(X, y, sample_names):
	'''
	I want training set to have equal number of each class.
	I want training to be .7 of total. Test gets .3
	'''

	numb_of_samps = len(y)

	numb_of_training = int(numb_of_samps * .7)
	#numb_of_training = 24

	numb_of_test = numb_of_samps - numb_of_training
	numb_of_class0_train = numb_of_training /2
	numb_of_class1_train = numb_of_class0_train

	#print numb_of_samps
	##print numb_of_training
	#print numb_of_test
	#print

	X_train =[]
	y_train =[]
	X_test =[]
	y_test =[]
	names_train = []
	names_test = []

	numb_of_class0_in_train_atm = 0
	numb_of_class1_in_train_atm = 0

	for i in range(len(y)):

		if y[i] == 0.0 and numb_of_class0_in_train_atm < numb_of_class0_train:
			X_train.append(X[i])
			y_train.append(y[i])
			numb_of_class0_in_train_atm += 1
			names_train.append(sample_names[i])
		elif y[i] == 1.0 and numb_of_class1_in_train_atm < numb_of_class1_train:
			X_train.append(X[i])
			y_train.append(y[i])
			numb_of_class1_in_train_atm += 1
			names_train.append(sample_names[i])
		else:
			X_test.append(X[i])
			y_test.append(y[i])
			names_test.append(sample_names[i])


	# print len(X_train)
	# print len(y_train)
	# print y_train
	# print len(X_test)
	# print len(y_test)
	# print y_test



	return np.array(X_train), np.array(y_train), np.array(X_test), np.array(y_test), names_train, names_test

