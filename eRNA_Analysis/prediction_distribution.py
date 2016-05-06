	
import numpy as np

#to plot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt



def prediction_distribution_plot(train_class0, test_class0, test_class1, train_class1, output_file_name, dot_or_name):


	if dot_or_name == 'dot':

		plt.figure(20)

		for i in train_class0:
				plt.scatter(1, i, c='blue')
		for i in test_class0:
				plt.scatter(2, i, c='lightblue')
		for i in test_class1:
				plt.scatter(3, i, c='pink')
		for i in train_class1:
				plt.scatter(4, i, c='red')
	

		plt.ylabel('Prediction Value')
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

		plt.savefig(output_file_name)




	if dot_or_name == 'blind_dot':

		plt.figure(20)

		for i in train_class0:
				plt.scatter(1, i, c='blue')
		for i in test_class0:
				plt.scatter(2, i, c='grey')
		for i in train_class1:
				plt.scatter(3, i, c='red')
	

		plt.ylabel('Prediction Value')
		plt.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off')

		plt.grid(b=True, axis='y', which='major', color='grey', linestyle='--')

		#plt.annotate('LG', (0.3,0.03), xycoords='figure fraction')
		#plt.annotate('HG', (0.7,0.03), xycoords='figure fraction')
		plt.annotate('LG Training', (0.2,0.07), xycoords='figure fraction')
		plt.annotate('Test', (0.49,0.07), xycoords='figure fraction')
		#plt.annotate('Test', (0.61,0.07), xycoords='figure fraction')
		plt.annotate('HG Training', (0.72,0.07), xycoords='figure fraction')

		#plt.xlim([-0.01, 1.01])
		plt.ylim([-0.01, 1.01])
		# x = np.arange(0, 4, 0.01)
		# plt.plot(x, [0.5]*len(x), '--', c='grey')

		plt.savefig(output_file_name)




	if dot_or_name == 'blind_name':

		plt.figure(20)

		for i in train_class0:
				plt.scatter(1, i, c='blue')
		for i in range(len(test_class0)):
				#plt.scatter(2, i, c='grey')
				plt.annotate(test_class1[i], xy=(2, test_class0[i]), size=4.9, bbox=dict(boxstyle="round", fc="0.8"), color='black')
		for i in train_class1:
				plt.scatter(3, i, c='red')
	

		plt.ylabel('Prediction Value')
		plt.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off')

		plt.grid(b=True, axis='y', which='major', color='grey', linestyle='--')

		#plt.annotate('LG', (0.3,0.03), xycoords='figure fraction')
		#plt.annotate('HG', (0.7,0.03), xycoords='figure fraction')
		plt.annotate('LG Training', (0.2,0.07), xycoords='figure fraction')
		plt.annotate('Test', (0.49,0.07), xycoords='figure fraction')
		#plt.annotate('Test', (0.61,0.07), xycoords='figure fraction')
		plt.annotate('HG Training', (0.72,0.07), xycoords='figure fraction')

		#plt.xlim([-0.01, 1.01])
		plt.ylim([-0.01, 1.01])
		# x = np.arange(0, 4, 0.01)
		# plt.plot(x, [0.5]*len(x), '--', c='grey')

		plt.savefig(output_file_name)










	# 	x_pos = 1.005
	# 	for old_point in results:
	# 		#if they are close together
	# 		if np.abs(score - old_point) < 0.04:
	# 			x_pos += 0.003

	# 	if this_samp in set3_LG:
	# 		plt.annotate(this_samp, xy=(x_pos, score), size=4.4, bbox=dict(boxstyle="round", fc="0.8"),  color='blue')
	# 	elif this_samp in set3_HG:
	# 		plt.annotate(this_samp, xy=(x_pos, score), size=4.4, bbox=dict(boxstyle="round", fc="0.8"), color='red')
	# 	else:
	# 		plt.annotate(this_samp, xy=(x_pos, score), size=3.9, bbox=dict(boxstyle="round", fc="0.8"), color='yellow')

	# 	results.append(score)















	# plt.vlines(1,0,1) 

	# #this_samp = sample_names[test_index]

	# x_pos = 1.005
	# for old_point in results:
	# 	#if they are close together
	# 	if np.abs(score - old_point) < 0.04:
	# 		x_pos += 0.003

	# if this_samp in set3_LG:
	# 	plt.annotate(this_samp, xy=(x_pos, score), size=4.4, bbox=dict(boxstyle="round", fc="0.8"),  color='blue')
	# elif this_samp in set3_HG:
	# 	plt.annotate(this_samp, xy=(x_pos, score), size=4.4, bbox=dict(boxstyle="round", fc="0.8"), color='red')
	# else:
	# 	plt.annotate(this_samp, xy=(x_pos, score), size=3.9, bbox=dict(boxstyle="round", fc="0.8"), color='yellow')

	# results.append(score)



	# #y = numpy.ones(numpy.shape(results))   # Make all y values the same
	# #plt.plot(y,results,'_',ms = 20)  # Plot a line at each location specified in a


	# acc = accuracy(targets, outputs)
	
	# x = np.arange(1, 1.07, 0.01)
	# #draws line at decision boundary (0.5)
	# plt.plot(x, [0.5]*len(x), '--', c='grey')

	# plt.xticks([])
	# plt.xlim(1,1.07)
	# plt.ylim(0,1.02)
	# plt.title("L1 Logistic Grade Predictions\neRNAs")
	# plt.ylabel("Score")
	# plt.annotate('HG', xy=(0.95,0.515), xycoords='axes fraction', size='small', bbox=dict(boxstyle="round", fc="0.8"), color='red')
	# plt.annotate('LG', xy=(0.95,0.45), xycoords='axes fraction', size='small', bbox=dict(boxstyle="round", fc="0.8"), color='blue')
	# plt.annotate('Accuracy= %0.2f' % acc, xy=(0.81,0.03), xycoords='axes fraction', size='small', bbox=dict(boxstyle="round", fc="0.8"))

	# #plt.annotate('Number of Genes = ' + str(numbOfSigGenes), xy=(0.7,0.7.5), xycoords='axes fraction', size='small', bbox=dict(boxstyle="round", fc="0.8"), color='maroon')

	# #plt.annotate('Testing Set = set_2', xy=(0.7,0.7), xycoords='axes fraction', size='small', bbox=dict(boxstyle="round", fc="0.8"), color='maroon')
	# #plt.annotate('Training Set = set_1', xy=(0.7,0.65), xycoords='axes fraction', size='small', bbox=dict(boxstyle="round", fc="0.8"), color='maroon')



	# plt.savefig('eRNA_predictions_1.pdf')