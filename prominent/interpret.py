import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import shap
import pickle
import seaborn as sns
import argparse

def get_feature_importance():
    parser = argparse.ArgumentParser(description='Caculate the importance for each feature and pathway.')
    parser.add_argument('--input_shap', default='shap.pkl', help='SHAP output from training.')
    parser.add_argument('--output_feature', default='Feature_importance.csv', help='Output feature importance file path. csv file.')
    parser.add_argument('--output_pathway', default='Pathway_importance.csv', help='Output pathway importance file path. csv file.')
    args = parser.parse_args()

    input_file_shap = args.input_shap
    output_file_feat = args.output_feature
    output_file_path = args.output_pathway
    
    with open(input_file_shap, 'rb') as file:
        shap_dict = pickle.load(file)
    shaps = shap_dict["path"]
    shaps_feat = shap_dict["feat"]
    features = pd.read_csv("pathway_gobp.csv",index_col=0).columns.tolist()
    result_array = np.zeros_like(np.sum(np.abs(shaps[0][1]), axis=0)[32:])
    pos_arr = np.zeros_like(np.sum(np.abs(shaps[0][1]), axis=0)[32:])
    neg_arr = np.zeros_like(np.sum(np.abs(shaps[0][1]), axis=0)[32:])
    for i, shap in enumerate(shaps):
        feature_importance = np.sum(np.abs(shap[1]), axis=0)[32:]
        positive_sum = np.sum(shap[1] * (shap[1] > 0), axis=0)[32:]
        negative_sum = np.sum(shap[1] * (shap[1] < 0), axis=0)[32:]
        result_array += feature_importance
        pos_arr += positive_sum
        neg_arr += negative_sum
    fi_result = pd.DataFrame({"path_name":features,"importance":result_array, "pos":pos_arr,"neg": np.abs(neg_arr)})
    df_sorted = fi_result.sort_values(by='importance', ascending=False)
    df_sorted.to_csv(output_file_path, index=False)

    features = pd.read_csv("gene.average.beta.by.intensity.csv",index_col=0).index.tolist()
    result_array = np.zeros_like(np.sum(np.abs(shaps_feat[0][1]), axis=0))
    pos_arr = np.zeros_like(np.sum(np.abs(shaps_feat[0][1]), axis=0))
    neg_arr = np.zeros_like(np.sum(np.abs(shaps_feat[0][1]), axis=0))
    for i, shap in enumerate(shaps_feat):
        feature_importance = np.sum(np.abs(shap[1]), axis=0)
        positive_sum = np.sum(shap[1] * (shap[1] > 0), axis=0)
        negative_sum = np.sum(shap[1] * (shap[1] < 0), axis=0)
        result_array += feature_importance
        pos_arr += positive_sum
        neg_arr += negative_sum
    fi_result = pd.DataFrame({"path_name":features,"importance":result_array, "pos":pos_arr,"neg": np.abs(neg_arr)})
    df_sorted = fi_result.sort_values(by='importance', ascending=False)
    df_sorted.to_csv(output_file_feat,index=False)
    
if __name__ == '__main__':
    get_feature_importance()