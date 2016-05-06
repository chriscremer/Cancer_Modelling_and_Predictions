

import numpy as np

from sklearn.metrics import roc_curve, auc, accuracy_score

import matplotlib.mlab as mlab
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def make_roc(y_valid, output, title, output_file_name):

	mean_fpr = np.linspace(0, 1, 100)
	
	#output = takata_fit(X_train, y_train, X_valid, gene_set, gene_names)
	fpr, tpr, thresholds = roc_curve(y_valid, output)
	roc_auc = auc(fpr, tpr)
	#print fpr
	#print tpr
	plt.plot(fpr, tpr, label='ROC curve (area = %0.2f)' % roc_auc)
	#roc_auc = auc(fpr, tpr)	
	#takata_auc.append(roc_auc)
	#mean_tpr_takata += interp(mean_fpr, fpr, tpr)
	#mean_tpr_takata[0] = 0.0

	# mean_tpr_takata /= len(takata_auc)
	# mean_tpr_takata[-1] = 1.0
	# mean_auc_takata = auc(mean_fpr, mean_tpr_takata)
	# plt.plot(mean_fpr, mean_tpr_takata, lw=1, label='Takata/Kato Gene Set (area = %0.2f)' % (mean_auc_takata))

	plt.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='Luck')
	#plt.plot(mean_fpr, mean_tpr_gnb, 'k--', label='Mean ROC (area = %0.2f)' % mean_auc, lw=2)

	plt.xlim([-0.01, 1.01])
	plt.ylim([-0.01, 1.01])
	plt.title(title)
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	#plt.title('Receiver operating characteristic')
	plt.legend(loc="lower right", prop={'size':8})
	plt.savefig(output_file_name)
	print 'Plot complete'




if __name__ == "__main__":

	y_valid = [0,0,0,1,1,1]
	output = [0,0,1,0,1,1]

	make_roc(y_valid, output)