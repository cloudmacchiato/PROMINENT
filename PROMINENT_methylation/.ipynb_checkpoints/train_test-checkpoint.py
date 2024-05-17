import os, torch
import numpy as np
import pandas as pd
from torch import nn
import torch.nn.functional as F
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt
import argparse
from PROMINENT_methylation.train import *

def preprocessing1(path_fn, gene_fn, label_fn):
    ### pathway datasets
    if os.path.exists(path_fn):
        pathway = pd.read_csv(path_fn, header=0)       
        print(">> Pathway Data :",path_fn)
        pathway_info = pathway.iloc[:,1:]
        pathway_info = pathway_info.values
        pathway_info = np.transpose (pathway_info)
        pathway_info = torch.FloatTensor(pathway_info)
        print("pathway matrix shape : ",pathway_info.shape)
        print("num_pathway : ",pathway_info.shape[0])
    else:
        pathway_info = None
        
    ### methylation datasets
    print(">> Methylation Gene-level Data:",gene_fn)
    data = pd.read_csv(gene_fn, header=0)

    expression = data.iloc[:,1:]
    features = data.iloc[:,0].tolist()
    expression = expression.values
    expression = np.transpose(expression)

    scaler = MinMaxScaler()
    scaler = scaler.fit(expression)
    expression = scaler.transform(expression)

    sample_dim = expression.shape[0]
    input_dim = expression.shape[1]

    #print dimension of sample and number of genes
    print("sample_dim : ",sample_dim)
    print("input_size (number of genes): ",input_dim)
    
    pheno = pd.read_csv(label_fn, index_col=0)
    status = np.array(pheno['label']).reshape(-1,1)

    patient = list(data.iloc[:,1:].columns.values.tolist()) 
    print("patient list : ",patient[0:6])
    print("feature list : ",features[0:6])
    
    return pathway_info, expression, status, features


def train():
    parser = argparse.ArgumentParser(description='5 fold CV')
    parser.add_argument('--input_csv', default='gene.average.beta.by.intensity.csv', help='Input csv file of X_i. Rows are genes, columns are samples.')
    parser.add_argument('--input_label', default='label.csv', help='Input label file path.')
    parser.add_argument('--input_path', default='pathway_gobp.csv', help='Pathway information from dataprep.')
    parser.add_argument('--output', default='pred.pkl', help='Prediction output file path. pkl file.')
    parser.add_argument('--output_shap', default='shap.pkl', help='Model interpretation output file path. pkl file.')
    parser.add_argument('--mlp', action='store_true', help='DNN only model.')

    args = parser.parse_args()

    input_file_csv = args.input_csv
    input_file_label = args.input_label
    input_file_path = args.input_path
    output_file = args.output
    output_file_shap = args.output_shap
    mlp = args.mlp
    
    pathway_info, x_data, status, features = preprocessing1(input_file_path, input_file_csv,input_file_label)
    
    trainArgs = {}
    trainArgs['x_data'] = x_data
    trainArgs['y_data'] = status
    trainArgs['pathway_info'] = pathway_info
    trainArgs['features'] = features
    trainArgs['num_fc_list'] = [32]
    trainArgs['lr_list'] = [0.001]
    trainArgs['device'] = '0'
    trainArgs['seed'] = 0
    trainArgs['filename'] = output_file
    trainArgs['filename2'] = output_file_shap
    
    if mlp:
        print("Run PROMINENT_DNN.")
        train = train_kfold_mlp(trainArgs)
    else:
        train = train_kfold2(trainArgs)
        
    result = train.kfold()
    
if __name__ == '__main__':
    train()
    
