

'''
This file actually plots the samples unlike the other that plots the scores of the samples.
The dimensions are the 2 eRNAs.
'''


#code for loading dataset
from load_eRNA_dataset import return_dataset_matrix


from sklearn import preprocessing
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





if __name__ == "__main__":

	dataset, sample_names, gene_names = return_dataset_matrix('/data1/morrislab/ccremer/eRNA/bc_eRNA_1k_noOverlap_41Samples_09.txt')

	labels = get_labels(sample_names)


	X = dataset
	y = labels


	#onlye use the two good eRNAs
	print gene_names[3836]
	print gene_names[4046]
	
	
	print X.shape
	X = X.T
	X = X[[3836, 4046]]
	X = X.T
	print X.shape
	
	
	plt.figure(1)

	for i in range(len(X)):
		if y[i] == 0.0:
			lg = plt.scatter(X[i][0], X[i][1], c='blue')
		else:
			hg = plt.scatter(X[i][0], X[i][1], c='red')


	plt.title("LG/HG Scatter Plot")
	plt.xlabel("eRNA6000 Expression")
	plt.ylabel("eRNA52064 Expression")
	


	plt.legend((lg, hg), 
				('Low Grade', 'High Grade'), 
				scatterpoints=1, 
				loc='lower right',  
				fontsize=8)

	plt.savefig('eRNA_scatter_plot.pdf')


