

import numpy as np

from sklearn.metrics import roc_curve, auc, accuracy_score

import matplotlib.mlab as mlab
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def make_roc(y_valid, output, output_file_name):

	

	



	# plt.figure(1)
	# plt.gca().set_position((.1, .3, .8, .6)) # to make a bit of room for extra text



	# plt.plot(fpr, tpr, label='AUC= %0.2f' % roc_auc)
	# plt.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6))
	# #plt.plot(mean_fpr, mean_tpr_gnb, 'k--', label='Mean ROC (area = %0.2f)' % mean_auc, lw=2)
	# plt.xlim([-0.01, 1.01])
	# plt.ylim([-0.01, 1.01])
	# #plt.title(title)
	# plt.xlabel('False Positive Rate')
	# plt.ylabel('True Positive Rate')
	# #plt.title('Receiver operating characteristic')
	# plt.legend(loc="lower right", prop={'size':8})

	# plt.figtext(.02, .02, "Leave-one-out cross-validation ROC curve of Logistic Regression with L1 regularization (AUC=0.97)", size='small')

	# plt.savefig(output_file_name)
	








	#print fpr
	#print tpr
	
	#roc_auc = auc(fpr, tpr)	
	#takata_auc.append(roc_auc)
	#mean_tpr_takata += interp(mean_fpr, fpr, tpr)
	#mean_tpr_takata[0] = 0.0

	# mean_tpr_takata /= len(takata_auc)
	# mean_tpr_takata[-1] = 1.0
	# mean_auc_takata = auc(mean_fpr, mean_tpr_takata)
	# plt.plot(mean_fpr, mean_tpr_takata, lw=1, label='Takata/Kato Gene Set (area = %0.2f)' % (mean_auc_takata))





	# fig = plt.figure(1)
	# ax1 = fig.add_axes((.1,.4,.8,.5))
	# ax1.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6))
	# #plt.plot(mean_fpr, mean_tpr_gnb, 'k--', label='Mean ROC (area = %0.2f)' % mean_auc, lw=2)
	# plt.xlim([-0.01, 1.01])
	# plt.ylim([-0.01, 1.01])
	# #plt.title(title)
	# plt.xlabel('False Positive Rate')
	# plt.ylabel('True Positive Rate')
	# #plt.title('Receiver operating characteristic')
	# plt.legend(loc="lower right", prop={'size':8})

	# fig.text(.1,.1,title)

	# fig.savefig(output_file_name)

	plt.figure(23)

	mean_fpr = np.linspace(0, 1, 100)
	fpr, tpr, thresholds = roc_curve(y_valid, output)
	roc_auc = auc(fpr, tpr)

	plt.plot(fpr, tpr, label='AUC= %0.2f' % roc_auc)
	plt.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6))
	#plt.plot(mean_fpr, mean_tpr_gnb, 'k--', label='Mean ROC (area = %0.2f)' % mean_auc, lw=2)
	plt.xlim([-0.01, 1.01])
	plt.ylim([-0.01, 1.01])
	#plt.title(title)
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	#plt.title('Receiver operating characteristic')
	plt.legend(loc="lower right", prop={'size':8})


	plt.savefig(output_file_name)

	print 'ROC Plot complete'




if __name__ == "__main__":

	y_valid = [0,0,0,1,1,1]
	output = [0,0,1,0,1,1]

	make_roc(y_valid, output, 'title', 'name.pdf')