


import numpy as np
import math
from scipy.stats import ttest_ind

def p_value(dataset, labels, feature_names):
	'''
	Return the significance of the means of the two classes for each feature
	'''

	#split the dataset into the two classes
	class_0_dataset = []
	class_1_dataset = []
	for i in range(len(dataset)):

		if labels[i] == 0.0:
			class_0_dataset.append(dataset[i])
		else:
			class_1_dataset.append(dataset[i])

	#make features the rows
	class_0_dataset = np.array(class_0_dataset).T
	class_1_dataset = np.array(class_1_dataset).T

	#values needed for the calculation
	numb_of_feats = len(dataset[0])
	numb_0_samps = len(class_0_dataset.T)
	numb_1_samps = len(class_1_dataset.T)


	#level of significance
	alpha = 0.001
	#Bonferroni Method (correction for multiple testing problem)
	corrected_alpha = alpha / numb_of_feats

	#calculate p-vlaues of each feature
	sig_feat_names = []
	sig_feat_p_value = []
	for i in range(len(dataset.T)):

		t_statistic, p_value = ttest_ind(class_0_dataset[i], class_1_dataset[i], equal_var=False)

		if p_value < corrected_alpha:
			sig_feat_p_value.append(p_value)
			sig_feat_names.append(feature_names[i])


		'''
		mean0 = np.mean(class_0_dataset[i])
		mean1 = np.mean(class_1_dataset[i])

		std0 = np.std(class_0_dataset[i])
		std1 = np.std(class_1_dataset[i])

		#variance divided by numb of samples
		v_n0 = (std0**2) / numb_0_samps
		v_n1 = (std1**2) / numb_1_samps

		#standard error
		SE = math.sqrt(v_n0 + v_n1)

		#degrees of freedom
		df = int(((v_n0 + v_n1)**2)  / (((v_n0**2) /(numb_0_samps-1)) + ((v_n1**2) /(numb_1_samps-1))))
		#df = min([numb_0_samps, numb_1_samps])
		#print df


		t_score = (mean0 - mean1) / SE

		p_value = x + y
		'''

	arg_sort = np.argsort(sig_feat_p_value)
	sorted_pvalues = sorted(sig_feat_p_value)
	for i in range(len(sig_feat_names)):
		print sig_feat_names[arg_sort[i]] + '  ' + str(sorted_pvalues[i])

	#print corrected_alpha





		




#if __name__ == "__main__":


