
import random
import numpy as np





def make_simulated(numb_of_samples, numb_of_features, numb_of_sig_features, numb_of_class0):

	data = [[] for i in xrange(numb_of_samples)]

	for i in range(numb_of_features - numb_of_sig_features):

		#mean = random.gauss(5,5)
		#std = random.gauss(2,2)

		for j in range(numb_of_samples):

			#data[j].append(random.gauss(mean, std))
			data[j].append(random.uniform(-3, 15))


	for i in range(numb_of_sig_features):

		class1_mean = random.gauss(5,5)
		class1_std = random.gauss(2,2)
		class2_mean = random.gauss(5,5)
		class2_std = random.gauss(2,2)

		for j in range(numb_of_samples):

			#if samp is even its class 0 else, class 1
			if j % 2 == 0:

				data[j].append(random.gauss(class1_mean, class1_std))

			else:
				data[j].append(random.gauss(class2_mean, class2_std))


	labels = []
	'''
	for i in range(numb_of_samples):
		if i % 2 == 0:
			labels.append(0.0)
		else:
			labels.append(1.0)
	'''
	for i in range(numb_of_class0):
		labels.append(0.0)
	for i in range(numb_of_samples - numb_of_class0):
		labels.append(1.0)


	return np.array(data), np.array(labels)





#a = make_simulated(10, 5, 2)
#print a[0]
#print a[1]