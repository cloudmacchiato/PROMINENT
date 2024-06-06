# PROMINENT
## Installation
`conda create -n prominent python=3.8`\
`conda activate prominent`\
`conda install pytorch torchvision`\
`pip install PROMINENT-methylation`

## Running PROMINENT
For helps, run `functions -h`.\
Data preparation for model input.\
`PROMINENT-data_prepare --input_csv <filename> --input_gmt <filename> --output <filename>`\
`PROMINENT-train_test_cv --input_csv <filename> --input_label <filename> --input_path <filename> --output <filename> --output_shap <filename> --mlp`\
    `--mlp` is optional for running PROMINENT_DNN\
`PROMINENT-scores --input_pkl <filename> --output <filename>`\
`PROMINENT-model_interpret --input_shap <filename> --output_feature <filename> --output_pathway <filename>`

