
import numpy as np
from sklearn.cross_validation import StratifiedKFold

#sklearn models
from sklearn import linear_model
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn import tree
from sklearn.neighbors import KNeighborsClassifier
from sklearn import lda
from sklearn.ensemble import AdaBoostClassifier


def accuracy(target, output, sample_names2, print_yes_no):
	'''Model accuracy'''

	results = []
	for i in range(len(output)):
		if output[i] > 0.5 and target[i] > 0.5:
			results.append(1.0)
		elif output[i] < 0.5 and target[i] < 0.5:
			results.append(1.0)
		else:
			results.append(0.0)
			if print_yes_no:
				print sample_names2[i]

	return np.mean(results)


def two_fold_selection(X_train, y_train, X_test, y_test, sample_names2):
	'''hyperparameter selection and test models'''

	gaussian_naive_bayes_acc = []
	mean_tpr_gnb = 0.0

	logistic_reg_L1_acc = []
	mean_tpr_lr1 = 0.0

	logistic_reg_L2_acc = []
	mean_tpr_lr2 = 0.0

	logistic_reg_elasticnet_acc = []
	mean_tpr_lrel = 0.0

	knn_acc = []
	mean_tpr_knn = 0.0

	random_forest_acc = []
	mean_tpr_rf = 0.0

	svc_acc = []
	mean_tpr_svc = 0.0

	lda_acc = []
	mean_tpr_lda = 0.0

	ada_acc = []
	mean_tpr_ada = 0.0

	mean_fpr = np.linspace(0, 1, 100)


	#hyperparameter selection
	cv_inner = StratifiedKFold(y_train, n_folds=2)
	for j, (train_index2, valid_index2) in enumerate(cv_inner):
		#print 'Selecting Good Hyperparameters'

		half1, half2 = X_train[train_index2], X_train[valid_index2]
		y1, y2 = y_train[train_index2], y_train[valid_index2]

		#print '1'
		#print 'len y1 ' + str(len(y1))
		#print 'len y2 ' + str(len(y2))

		'''
		best_log_reg_L1_auc = 0
		best_log_reg_L1_hyper = 0
		for hyper in [0.001, 0.01, 0.05, 0.1, 1.0, 2.0]:
			auc_sum = 0
			for k in [0,1]:
				if k == 0:
					first_half = half1
					first_y = y1
					other_half = half2
					other_y = y2
				else:
					first_half = half2
					first_y = y2
					other_half = half1
					other_y= y1				

				clf = linear_model.LogisticRegression(penalty='l1', C=hyper)
				#clf = linear_model.Lasso(alpha=hyper)
				output = clf.fit(first_half, first_y).predict(other_half)
				acc = accuracy(other_y, output, sample_names2)
				#print hyper
				#print output
				#print other_y
				#print acc

				#fpr, tpr, thresholds = roc_curve(other_y, output[:, 0])
				#roc_auc = auc(fpr, tpr)
				auc_sum += acc
			#print 'L1 hyper ' + str(hyper)
			#print 'auc_sum ' + str(auc_sum)
			if auc_sum > best_log_reg_L1_auc:
				best_log_reg_L1_auc = auc_sum
				best_log_reg_L1_hyper = hyper
		#print 'Best L1 hyper ' + str(best_log_reg_L1_hyper)
		'''




		#print '2'


		best_log_reg_L2_auc = 0
		best_log_reg_L2_hyper = 0
		for hyper in [0.0001, 0.001, 0.01, 0.1, 1.0, 10]:
			auc_sum = 0
			for k in [0,1]:
				if k == 0:
					first_half = half1
					first_y = y1
					other_half = half2
					other_y = y2
				else:
					first_half = half2
					first_y = y2
					other_half = half1
					other_y= y1				

				clf = linear_model.LogisticRegression(penalty='l2', C=hyper)
				output = clf.fit(first_half, first_y).predict(other_half)
				acc = accuracy(other_y, output, sample_names2, False)

				#fpr, tpr, thresholds = roc_curve(other_y, output[:, 0])
				#roc_auc = auc(fpr, tpr)
				auc_sum += acc
			if auc_sum > best_log_reg_L2_auc:
				best_log_reg_L2_auc = auc_sum
				best_log_reg_L2_hyper = hyper
		#print 'Best L2 hyper ' + str(best_log_reg_L2_hyper)


		'''
		print '3'

		best_log_reg_el_auc = 0
		best_log_reg_el_hyper = (0,0)
		for l1_ratio in [0.01, 0.25, 0.5, 0.75, 1.0]:
			for alpha in [0.01, 0.1, 1.0, 10]:
				auc_sum = 0
				for k in [0,1]:
					if k == 0:
						first_half = half1
						first_y = y1
						other_half = half2
						other_y = y2
					else:
						first_half = half2
						first_y = y2
						other_half = half1
						other_y= y1				

					clf = linear_model.ElasticNet(l1_ratio=l1_ratio, alpha=alpha)
					output = clf.fit(first_half, first_y).predict(other_half)
					acc = accuracy(y_test, output)

					#fpr, tpr, thresholds = roc_curve(other_y, output[:, 0])
					#roc_auc = auc(fpr, tpr)
					auc_sum += acc
				if auc_sum > best_log_reg_el_auc:
					best_log_reg_el_auc = auc_sum
					best_log_reg_el_hyper = (l1_ratio, alpha)
		#print 'Best EL hyper ' + str(best_log_reg_el_hyper)
		'''



		#print '4'

		'''
		best_knn_auc = 0
		best_knn_hyper = 0
		for hyper in range(1,5,2):
			auc_sum = 0
			for k in [0,1]:
				if k == 0:
					first_half = half1
					first_y = y1
					other_half = half2
					other_y = y2
				else:
					first_half = half2
					first_y = y2
					other_half = half1
					other_y= y1				

				clf = KNeighborsClassifier(n_neighbors=hyper)
				output = clf.fit(first_half, first_y).predict(other_half)
				acc = accuracy(other_y, output, sample_names2)

				#fpr, tpr, thresholds = roc_curve(other_y, output[:, 0])
				#roc_auc = auc(fpr, tpr)
				auc_sum += acc
			if auc_sum > best_knn_auc:
				best_knn_auc = auc_sum
				best_knn_hyper = hyper
		#print 'Best k hyper ' + str(best_knn_hyper)



		#print '5'

		best_rf_auc = 0
		best_rf_hyper = 0
		for hyper in range(1,152,50):
			auc_sum = 0
			for k in [0,1]:
				if k == 0:
					first_half = half1
					first_y = y1
					other_half = half2
					other_y = y2
				else:
					first_half = half2
					first_y = y2
					other_half = half1
					other_y= y1				

				clf = RandomForestClassifier(n_estimators=hyper)
				output = clf.fit(first_half, first_y).predict(other_half)
				acc = accuracy(other_y, output, sample_names2)

				#fpr, tpr, thresholds = roc_curve(other_y, output[:, 0])
				#roc_auc = auc(fpr, tpr)
				auc_sum += acc
			if auc_sum > best_rf_auc:
				best_rf_auc = auc_sum
				best_rf_hyper = hyper
		#print 'Best k hyper ' + str(best_knn_hyper)


		#print '6'

		best_svc_auc = 0
		best_svc_hyper = 0
		for hyper in [0.0001, 0.001, 0.01, 0.1, 1.0, 10]:
			auc_sum = 0
			for k in [0,1]:
				if k == 0:
					first_half = half1
					first_y = y1
					other_half = half2
					other_y = y2
				else:
					first_half = half2
					first_y = y2
					other_half = half1
					other_y= y1				

				clf = SVC(C=hyper, probability=True)
				output = clf.fit(first_half, first_y).predict(other_half)
				acc = accuracy(other_y, output, sample_names2)

				#fpr, tpr, thresholds = roc_curve(other_y, output[:, 0])
				#roc_auc = auc(fpr, tpr)
				auc_sum += acc
			if auc_sum > best_svc_auc:
				best_svc_auc = auc_sum
				best_svc_hyper = hyper
		#print 'Best L1 hyper ' + str(best_log_reg_L1_hyper)



		#print '7'

		best_ada_auc = 0
		best_ada_hyper = 0
		for hyper in [2, 10, 50, 100, 200]:
			auc_sum = 0
			for k in [0,1]:
				if k == 0:
					first_half = half1
					first_y = y1
					other_half = half2
					other_y = y2
				else:
					first_half = half2
					first_y = y2
					other_half = half1
					other_y= y1				

				clf = AdaBoostClassifier(n_estimators=hyper)
				output = clf.fit(first_half, first_y).predict(other_half)
				acc = accuracy(other_y, output, sample_names2)

				#fpr, tpr, thresholds = roc_curve(other_y, output[:, 0])
				#roc_auc = auc(fpr, tpr)
				auc_sum += acc
			if auc_sum > best_ada_auc:
				best_ada_auc = auc_sum
				best_ada_hyper = hyper
		#print 'Best L1 hyper ' + str(best_log_reg_L1_hyper)
		'''

		


		#since only 2 fold, I run on both halves already to get average
		#hyperparameter selection complete
		break


	#Test Models with their best hyperparameters
	'''
	clf = GaussianNB()
	output = clf.fit(X_train, y_train).predict(X_test)
	#AUC
	#fpr, tpr, thresholds = roc_curve(y_test, output)
	#roc_auc = auc(fpr, tpr)
	#mean_tpr_gnb += interp(mean_fpr, fpr, tpr)
	#mean_tpr_gnb[0] = 0.0
	#Accuracy
	acc = accuracy(y_test, output, sample_names2)
	gaussian_naive_bayes_acc.append(acc)


	clf = lda.LDA()
	output = clf.fit(X_train, y_train).predict(X_test)
	#AUC
	#fpr, tpr, thresholds = roc_curve(y_test, output)
	#roc_auc = auc(fpr, tpr)	
	#mean_tpr_lda += interp(mean_fpr, fpr, tpr)
	#mean_tpr_lda[0] = 0.0
	#Accuracy
	acc = accuracy(y_test, output, sample_names2)
	lda_acc.append(acc)

	clf = linear_model.LogisticRegression(penalty='l1', C=best_log_reg_L1_hyper)
	#print 'best L1 ' + str(best_log_reg_L1_hyper)
	#clf = linear_model.Lasso(alpha=best_log_reg_L1_hyper)
	output = clf.fit(X_train, y_train).predict(X_test)
	#AUC
	#fpr, tpr, thresholds = roc_curve(y_test, output)
	#roc_auc = auc(fpr, tpr)	
	#mean_tpr_lr1 += interp(mean_fpr, fpr, tpr)
	#mean_tpr_lr1[0] = 0.0
	#Accuracy
	acc = accuracy(y_test, output, sample_names2)
	logistic_reg_L1_acc.append(acc)
	'''

	clf = linear_model.LogisticRegression(penalty='l2', C=best_log_reg_L2_hyper)
	output = clf.fit(X_train, y_train).predict(X_test)
	#AUC
	#fpr, tpr, thresholds = roc_curve(y_test, output)
	#roc_auc = auc(fpr, tpr)	
	#mean_tpr_lr2 += interp(mean_fpr, fpr, tpr)
	#mean_tpr_lr2[0] = 0.0
	#Accuracy
	acc = accuracy(y_test, output, sample_names2, True)
	logistic_reg_L2_acc.append(acc)

	'''
	clf = linear_model.ElasticNet(l1_ratio=best_log_reg_el_hyper[0], alpha=best_log_reg_el_hyper[1])
	output = clf.fit(X_train, y_train).predict(X_test)
	#AUC
	#fpr, tpr, thresholds = roc_curve(y_test, output)
	#roc_auc = auc(fpr, tpr)	
	#mean_tpr_lrel += interp(mean_fpr, fpr, tpr)
	#mean_tpr_lrel[0] = 0.0
	#Accuracy
	acc = accuracy(y_test, output)
	logistic_reg_elasticnet_acc.append(acc)
	'''

	'''
	clf = KNeighborsClassifier(n_neighbors=best_knn_hyper)
	output = clf.fit(X_train, y_train).predict(X_test)
	#AUC
	#fpr, tpr, thresholds = roc_curve(y_test, output)
	#roc_auc = auc(fpr, tpr)	
	#mean_tpr_knn += interp(mean_fpr, fpr, tpr)
	#mean_tpr_knn[0] = 0.0
	#Accuracy
	acc = accuracy(y_test, output, sample_names2)
	knn_acc.append(acc)


	clf = RandomForestClassifier(n_estimators=best_rf_hyper)
	output = clf.fit(X_train, y_train).predict(X_test)
	#AUC
	#fpr, tpr, thresholds = roc_curve(y_test, output)
	#roc_auc = auc(fpr, tpr)	
	#mean_tpr_rf += interp(mean_fpr, fpr, tpr)
	#mean_tpr_rf[0] = 0.0
	#Accuracy
	acc = accuracy(y_test, output, sample_names2)
	random_forest_acc.append(acc)


	clf = SVC(C=best_svc_hyper, probability=True)
	output = clf.fit(X_train, y_train).predict(X_test)
	#AUC
	#fpr, tpr, thresholds = roc_curve(y_test, output)
	#roc_auc = auc(fpr, tpr)	
	#mean_tpr_svc += interp(mean_fpr, fpr, tpr)
	#mean_tpr_svc[0] = 0.0
	#Accuracy
	acc = accuracy(y_test, output, sample_names2)
	svc_acc.append(acc)


	clf = AdaBoostClassifier(n_estimators=best_ada_hyper)
	output = clf.fit(X_train, y_train).predict(X_test)
	#AUC
	#fpr, tpr, thresholds = roc_curve(y_test, output)
	#roc_auc = auc(fpr, tpr)	
	#mean_tpr_ada += interp(mean_fpr, fpr, tpr)
	#mean_tpr_ada[0] = 0.0
	#Accuracy
	acc = accuracy(y_test, output, sample_names2)
	ada_acc.append(acc)
	'''

	return (['GNB', 'LDA', 'LR1', 'LR2', 'Elastic', 'KNN', 'RandForest', 'SVC', 'AdaBoost'], 

			[np.mean(gaussian_naive_bayes_acc), 
			np.mean(lda_acc),
			np.mean(logistic_reg_L1_acc), 
			np.mean(logistic_reg_L2_acc), 
			np.mean(logistic_reg_elasticnet_acc), 
			np.mean(knn_acc), 
			np.mean(random_forest_acc), 
			np.mean(svc_acc),
			np.mean(ada_acc),
			])
