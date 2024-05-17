import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import shap
import pickle
import seaborn as sns
from sklearn.metrics import roc_curve, precision_recall_curve, auc
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.metrics import precision_recall_curve, matthews_corrcoef,accuracy_score
from scipy.stats import shapiro,mannwhitneyu,ttest_ind
import seaborn as sns
import argparse

def softmax(x):
    # Subtract the maximum value from each element for numerical stability
    x_exp = np.exp(x - np.max(x, axis=1, keepdims=True))
    # Calculate the softmax probabilities
    return x_exp / np.sum(x_exp, axis=1, keepdims=True)
def fscore(true, prob, beta):
    precision,recall,_ = precision_recall_curve(true,prob)
    f = (1+beta**2)*(precision*recall)/(beta**2*precision+recall)
    idx = np.nanargmax(f)
    pr = precision[idx] 
    rc = recall[idx] 
    f = f[idx]
    return pr, rc, f
def calculate_metrics(y_true, y_pred):
    y_true = np.array(y_true).astype(int)
    TP = np.sum((y_pred == 1) & (y_true == 1))
    TN = np.sum((y_pred == 0) & (y_true == 0))
    FP = np.sum((y_pred == 1) & (y_true == 0))
    FN = np.sum((y_pred == 0) & (y_true == 1))
    TPR = TP / (TP + FN)
    TNR = TN / (TN + FP)
    PPV = TP / (TP + FP)
    NPV = TN / (TN + FN)
    return TPR, TNR, PPV, NPV
def mcc_score(y_true, y_prob):
    mccs = []
    thresholds = np.arange(0,1,0.001)
    for threshold in thresholds:
        y_pred = np.greater(y_prob, threshold).astype(int)
        mcc = matthews_corrcoef(y_true, y_pred)
        mccs.append(mcc)
    mccs = np.array(mccs)
    max_mcc = mccs.max() 
    max_mcc_threshold =  thresholds[mccs.argmax()]
    return max_mcc, max_mcc_threshold

def get_scores():
    parser = argparse.ArgumentParser(description='Caculate classification scores.')
    parser.add_argument('--input_pkl', default='pred.pkl', help='Precition output file from training. Default: pred.pkl')
    parser.add_argument('--output', default='scores.csv', help='Output file path. csv file.')

    args = parser.parse_args()

    input_file_pkl = args.input_pkl
    output_file = args.output
    
    with open(input_file_pkl, 'rb') as file:
        test_preds = pickle.load(file)
    scores = []
    for i in range(0,5):
        y_true = np.array(test_preds[i][0]).astype(int)
        y_prob = np.array(test_preds[i][1])
        auc = roc_auc_score(y_true,y_prob)
        prauc = average_precision_score(y_true,y_prob)
        pr1, rc1, f1 = fscore(y_true,y_prob,1)
        pr2, rc2, f2 = fscore(y_true,y_prob,2)
        pr5, rc5, f5 = fscore(y_true,y_prob,5)
        pr10, rc10, f10 = fscore(y_true,y_prob,10)
        mcc, threshold = mcc_score(y_true, y_prob)
        y_pred = np.greater(y_prob, threshold).astype(int)
        acc = accuracy_score(y_true, y_pred)
        tpr, tnr, ppv, npv = calculate_metrics(y_true, y_pred)
        scores.append([auc,prauc, mcc,threshold, acc,tpr, tnr, ppv, npv,f1,f2,f5,f10])
    data = np.array(scores)
    row_names = ['Fold1', 'Fold2', 'Fold3', 'Fold4', 'Fold5']
    col_names = ["AUROC","PRAUC", "MCC", "Threshold(by MCC)", "Accuracy","TPR(recall)", "TNR(specificity)", "PPV(precision)", "NPV","Fscore(beta=1)", "Fscore(beta=2)", "Fscore(beta=5)","Fscore(beta=10)"]
    df = pd.DataFrame(data, index=row_names, columns=col_names)
    
    df.to_csv(output_file)
    
if __name__ == '__main__':
    get_scores()
    