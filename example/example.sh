#!/bin/sh
unzip example_data.zip

# To get helps, run -h or --help. For example:
PROMINENT-train_test --help
# Data preparation
PROMINENT-data_prepare --input_csv gene.average.beta.by.intensity.csv --input_gmt c5.go.bp.v2023.1.Hs.symbols.gmt --output pathway_gobp.csv 
#5-fold CV, PROMINENT_GOBP
PROMINENT-train_test_cv --input_csv gene.average.beta.by.intensity.csv --input_label label.csv --input_path pathway_gobp.csv --output pred.pkl --output_shap shap.pkl
# For PROMINENT_DNN
PROMINENT-train_test --input_csv gene.average.beta.by.intensity.csv --input_label label.csv --input_path pathway_gobp.csv --output pred.pkl --output_shap shap.pkl --mlp
# Classification evaluation
PROMINENT-scores --input_pkl pred.pkl --output scores.csv 
# Get feature and pathway importance
PROMINENT-model_interpret --input_shap shap.pkl --output_feature Feature_importance.csv --output_pathway Pathway_importance.csv
# For two independent datasets (--mlp if you want to use PROMINENT_DNN):
unzip example_two_datasets.zip
PROMINENT-model_train_test_independent --input_train_csv gene.average.beta.by.intensity.train.csv --input_test_csv gene.average.beta.by.intensity.test.csv --input_label_train label_train.csv --input_label_test label_test.csv --input_pathway pathway_gobp.csv
