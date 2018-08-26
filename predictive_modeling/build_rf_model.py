import os
import sys
import numpy as np
import random
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingRegressor
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.decomposition import PCA
from sklearn.neural_network import MLPClassifier
from sklearn.naive_bayes import BernoulliNB
from sklearn.naive_bayes import GaussianNB
from sklearn import metrics

positives = sys.argv[1]
negatives = sys.argv[2]
estimator_count = sys.argv[3]
max_depth = sys.argv[4]


def plot_ROC_curve(y_true, y_score, path, title):
	""" Plot and save the ROC curve by varying a threshold on y_score.
	Assume that the positive class is labeled truly as 1, and that if given a threshold value,
	we can classify the examples with: (y_score > threshold).astype(int)
	Return:
		The area under your curve (AUROC)
		The threshold value to use that would give you the highest F-score. Remember that each
		  point in an ROC curve has an associated F-score
	"""
	# print(y_score)
	def calc_f_score(thresh,y_true,y_score):
		# y_score = y_score.tolist()
		y_score_bin = y_score >= thresh

		# print(y_score_bin)
		TP = 0
		FP = 0
		FN = 0
		TN = 0
		for x in range(len(y_true)):
			if y_true[x] == True:
				if y_score_bin[x] == True:
					TP += 1
				elif y_score_bin[x] == False:
					FN += 1
			elif y_true[x] == False:
				if y_score_bin[x] == True:
					FP += 1
				elif y_score_bin[x] == False:
					TN += 1
		# print(TP,FP,FN,TN)
		TPR = TP/(TP + FN)
		FPR = 1 - TN/(TN + FP)
		f1 =  (2*TP)/((2*TP)+FP+FN)
		# if f1 > 0.95:
		# 	f1 = 0
			# f1 = 0
		return f1, TPR, FPR

	def find_max_ind(i_list):
		# print(i_list)
		cur_max = -0.1
		cur_index = -1
		for x in range(len(i_list)):
			if i_list[x] > cur_max:
				cur_max = i_list[x]
				cur_index = x
		return cur_index

	y_score_list = y_score.tolist()
	y_true_list = y_true.tolist()
	zipped_s = sorted(zip(y_score_list,y_true_list))
	# print(zipped_s)
	# sorted_ys = zipped
	thresholds = []
	f_scores = []
	TPs = []
	FPs = []
	current_thresh = -0.1

	for ex in zipped_s:
		if ex[0] > current_thresh:
			thresholds.append(ex[0])
			current_thresh = ex[0]
			f1,TP, FP = calc_f_score(current_thresh, y_true, y_score)
			TPs.append(TP)
			FPs.append(FP)
			f_scores.append(f1)

	auc = metrics.auc(FPs, TPs)
	print(f_scores)
	fig, ax = plt.subplots()
	plt.plot(FPs,TPs)
	plt.xlabel('FPR')
	plt.ylabel('TPR')
	plt.title(title + ", AUC="+str(auc)[0:6])

	plt.savefig(path)
	opt_ind = find_max_ind(f_scores)
	print(opt_ind)
	f_score_optimal_thr = thresholds[opt_ind]

	return auc, f_score_optimal_thr

X = []
Y = []

with open(positives, 'r') as pos:
    for line in pos:
        split_l = line.rstrip('\n').split('\t')
        values = split_l[3:]
        values_num = [float(x) if float(x) < 1000 else 1000 for x in values]
        X.append(values_num)
        Y.append(1)

with open(negatives, 'r') as neg:
    for line in neg:
        split_l = line.rstrip('\n').split('\t')
        values = split_l[3:]
        values_num = [float(x) if float(x) < 1000 else 1000 for x in values]
        X.append(values_num)
        Y.append(0)


inds = list(range(len(X)))
random.shuffle(inds)
X = [X[y] for y in inds]
Y = [Y[y] for y in inds]

train_cutoff = int((len(X) * 0.7)// 1)
valid_cutoff = int(train_cutoff + (len(X)*0.3*0.8)//1)
print(train_cutoff, valid_cutoff)
train_X = X[0:train_cutoff]
train_Y = Y[0:train_cutoff]
valid_X = X[train_cutoff:valid_cutoff]
# print(valid_X)
valid_Y = Y[train_cutoff:valid_cutoff]
test_X = X[valid_cutoff:]
test_Y = Y[valid_cutoff:]

model = GradientBoostingRegressor(n_estimators = int(estimator_count), max_depth=int(max_depth))
print(train_X)
print(train_Y)
model.fit(train_X, train_Y)

train_pred = model.predict(train_X)
# print(train_pred)
valid_pred = model.predict(valid_X)

plot_ROC_curve(np.array(train_Y), np.array(train_pred), 'training_roc_GBR'+estimator_count+"_"+max_depth+".png", 'ROC Curve, Training, Estimators='+estimator_count+", Max_Depth="+max_depth)

plot_ROC_curve(np.array(valid_Y), np.array(valid_pred), 'valid_roc_GBR'+estimator_count+"_"+max_depth+'.png', 'ROC Curve, Validation, Estimators='+estimator_count+", Max_Depth="+max_depth)

print('features: ')
print(model.feature_importances_)



