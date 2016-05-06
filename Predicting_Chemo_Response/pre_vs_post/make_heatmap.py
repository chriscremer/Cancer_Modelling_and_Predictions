

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np



def make_heatmap_func(matrix, x_ticks_names, plot_name):

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

	plt.xticks(np.arange(0,len(matrix[0]))+0.5, x_ticks_names, rotation='vertical', fontsize=5)

	#plt.xticks(np.arange(0,23)+0.5, x_ticks_names[:23], rotation='vertical', fontsize=5, color='blue')
	#plt.xticks(np.arange(23,len(matrix[0]))+0.5, x_ticks_names[23:], rotation='vertical', fontsize=5, color='red')
	
	plt.yticks(np.arange(0,len(matrix))+0.5, ['PER1'], fontsize=8)

	plt.tick_params(axis='both', which='both', left='off', right='off', bottom='off', top='off')
	
	plt.title('Pre vs Post Chemo\nTURBT vs RC', fontsize='small')

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



	plt.savefig(plot_name)





def make_heatmap_func_2(matrix, patient_ids, plot_name):

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

	plt.subplot(121)

	plt.pcolor(matrix, cmap=my_cmap, alpha=0.8, edgecolors='k', linewidths=0.5)

	plt.xlim([0,len(matrix[0])])

	plt.xticks(np.arange(0,len(matrix[0]))+0.5, ['TURBT', 'RC'], fontsize=8)

	#plt.xticks(np.arange(0,23)+0.5, x_ticks_names[:23], rotation='vertical', fontsize=5, color='blue')
	#plt.xticks(np.arange(23,len(matrix[0]))+0.5, x_ticks_names[23:], rotation='vertical', fontsize=5, color='red')
	
	plt.yticks(np.arange(0,len(matrix))+0.5, patient_ids, fontsize=8)

	plt.tick_params(axis='both', which='both', left='off', right='off', bottom='off', top='off')
	
	plt.title('PER1 Pre vs Post Chemo', fontsize='small')

	plt.colorbar(shrink=0.5)
	#plt.colorbar()


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



	plt.savefig(plot_name)