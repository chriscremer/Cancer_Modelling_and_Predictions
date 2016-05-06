

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np

from sklearn import preprocessing



def make_heatmap_func(matrix, title, x_ticks_names, y_ticks_names, output_file_name, weights=None, scale_the_features=False, sort_the_samples=False, labels=None):
	'''
	Given a matrix, turn it into a HEATMAP
	'''

	if sort_the_samples:


		LG_HG_sorted_indexes = np.argsort(labels)

		#for i in LG_HG_sorted_indexes:
		#	print str(LG_HG_sorted_indexes[i]) + '  ' + str(labels[i])

		matrix = matrix.T[LG_HG_sorted_indexes].T
		#labels = labels[LG_HG_sorted_indexes]
		new_sample_names = []
		for i in LG_HG_sorted_indexes:
			new_sample_names.append(x_ticks_names[i])

		x_ticks_names = new_sample_names



	#sort eRNAs

	weights_sorted_indexes = np.argsort(weights)

	#for i in LG_HG_sorted_indexes:
	#	print str(LG_HG_sorted_indexes[i]) + '  ' + str(labels[i])

	matrix = matrix[weights_sorted_indexes]
	#labels = labels[LG_HG_sorted_indexes]
	new_name_order = []
	for i in weights_sorted_indexes:
		new_name_order.append(y_ticks_names[i])

	y_ticks_names = new_name_order




	# #print X

	# #scale all expressions

	#make samps the rows
	matrix_T = matrix.T

	preprocessor = preprocessing.MinMaxScaler()
	preprocessor.fit(matrix_T)
	matrix_T = preprocessor.transform(matrix_T)
	matrix = matrix_T.T



	#HEATMAP

	#print matrix.shape
	#print len(matrix)
	#print len(matrix[1])
	#print len(x_ticks_names)
	#print matrix

	plt.figure(1)


	cdict = {'red' : ((0., 1., 1.), (0.45, 0., 0.), (1., 0., 0.)),
			'green': ((0., 0., 0.), (0.55, 0., 0.), (1., 1., 1.)), 
			'blue' : ((0., 0., 0.), (0.9, 0., 0.), (1., 0., 0.))
			}

	#cdict = {'red' : ((0., 0., 0.), (0.45, 0., 0.), (1., 0., 0.)),
	#		'green': ((0., 0., 0.), (0.5, 0.5, 0.5), (1., 1., 1.)), 
	#		'blue' : ((0., 0., 0.), (0.9, 0., 0.), (1., 0., 0.))
	#		}

	my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,6)

	#list_of_samps = ['resp', 'resp', 'resp', 'resp', 'resp', 'resp', 'resp', 'resp', 
	#				'non_resp',	'non_resp',	'non_resp',	'non_resp',	'non_resp',	'non_resp',	'non_resp',	'non_resp',	]
	#make heat map of significant genes/samples

	plt.subplot(211)

	plt.pcolor(matrix, cmap=my_cmap, alpha=0.8, edgecolors='k', linewidths=0.5)

	plt.xlim([0,len(matrix[0])])
	plt.ylim([0,len(matrix)])

	plt.xticks(np.arange(0,len(matrix[0]))+0.5, x_ticks_names, rotation='vertical', fontsize=5)

	#plt.xticks(np.arange(0,23)+0.5, x_ticks_names[:23], rotation='vertical', fontsize=5, color='blue')
	#plt.xticks(np.arange(23,len(matrix[0]))+0.5, x_ticks_names[23:], rotation='vertical', fontsize=5, color='red')
	
	plt.yticks(np.arange(0,len(matrix))+0.5, y_ticks_names, fontsize=5)

	plt.tick_params(axis='both', which='both', left='off', right='off', bottom='off', top='off')
	
	#plt.title(title, fontsize='small')

	#plt.colorbar(shrink=0.5)
	plt.colorbar()


	#plt.text(0.5, 0.5, 'matplotlib', horizontalalignment='center', verticalalignment='center')
	#plt.annotate( ' yo', xy=(2, -1))#, xycoords='axes fraction')



	#list_of_samps = ['resp', 'resp', 'resp', 'resp', 'resp', 'resp', 'resp',  
	#			'non_resp',	'non_resp',	'non_resp',	'non_resp',	'non_resp',	'non_resp', 'non_resp']
	#plot test expressions to the right of the training

	'''
	plt.subplot(122)
	plt.pcolor(test_both_matrix, cmap=my_cmap, alpha=0.8)
	plt.xticks(np.arange(0,len(both_matrix[0]))+0.5,test_samp_names, rotation='vertical', fontsize=5)
	#plt.yticks(np.arange(0,12)+0.5, top_gene_names, fontsize=8)
	plt.tick_params(axis='both', which='both', left='off', right='off', bottom='off', top='off', labelleft='off')
	plt.title('Expression Pattern of Distinguishing\nGenes in Test Data', fontsize='small')

	'''



	#plt.annotate('Responders', (0.24,0.014), xycoords='figure fraction', fontsize=7)
	#plt.annotate('Non-Responders', (0.43,0.014), xycoords='figure fraction', fontsize=7)
	#plt.annotate('Responders', (0.63,0.014), xycoords='figure fraction', fontsize=7)
	#plt.annotate('Non-Responders', (0.76,0.014), xycoords='figure fraction', fontsize=7)

	plt.tight_layout()



	plt.savefig(output_file_name)