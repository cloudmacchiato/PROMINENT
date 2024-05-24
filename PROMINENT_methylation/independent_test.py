import os, torch
import numpy as np
import pandas as pd
from torch import nn
import torch.nn.functional as F
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt
import argparse
from sklearn.metrics import roc_curve, precision_recall_curve, auc
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.metrics import precision_recall_curve, matthews_corrcoef,accuracy_score

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
def get_threshold(y_true_train, y_prob_train):
    mccs = []
    thresholds = np.arange(0,1,0.001)
    for threshold in thresholds:
        y_pred_train = np.greater(y_prob_train, threshold).astype(int)
        mcc = matthews_corrcoef(y_true_train, y_pred_train)
        mccs.append(mcc)
    mccs = np.array(mccs)
    max_mcc = mccs.max() 
    max_mcc_threshold =  thresholds[mccs.argmax()]
    return max_mcc_threshold    

def mcc_score(y_true, y_prob,threshold):
    y_pred = np.greater(y_prob, threshold).astype(int)
    mcc = matthews_corrcoef(y_true, y_pred) 
    return mcc, threshold

class train_test:
    def __init__(self, trainArgs):
   
        """
        args:
            x_data          : gene expression data
            y_data          : patient label
            pathway_info    : pathway matrics
            num_fc_list     : number of fully connected nodes 
            lr_list         : learning rate
            device          : GPU device
        Returns:
            AUC, Precision, Recall, F1
        """
        self.trainArgs = trainArgs
        self.seed_worker(trainArgs['seed'])
        os.environ["CUDA_VISIBLE_DEVICES"]=trainArgs['device']
        self.device = torch.device('cuda')
        directory = './CV_best_model'
        if not os.path.exists(directory):
            os.makedirs(directory)

    def kfold(self):
        trainArgs = self.trainArgs
        input_dim =expression.shape[1]
        pathway_info = trainArgs['pathway_info'].to(self.device)
        num_pathway = pathway_info.shape[0]
        num_fc_list = trainArgs['num_fc_list']
        lr_list = trainArgs['lr_list']
        random_seed = trainArgs['seed']
       
        kfold = StratifiedKFold(n_splits = 10, shuffle=True, random_state = random_seed)
        best_test_auc = 0
        best_fold = 0
        test_dataset1 = CustomDataset(expression1,status1) 
        test_loader1 = DataLoader(dataset = test_dataset1, batch_size = 64, shuffle = False)
        for fold, (train_index, test_index) in enumerate(kfold.split(expression, status)):   
            print('****************************************************************************')
            print('Fold {} / {}'.format(fold + 1 , kfold.get_n_splits()))
            print('****************************************************************************')
            x_train_ = expression[train_index]
            y_train_ = status[train_index] 
            x_test = expression[test_index]  
            y_test = status[test_index] 
            x_train, x_val, y_train, y_val = train_test_split(x_train_, y_train_, 
                                                              test_size=1/9, random_state = random_seed, stratify = y_train_)
            # smote = SMOTE(random_state=random_seed)
            # x_train_ = expression
            # y_train_ = status
            # x_test = expression1
            # y_test = status1
            # x_train, x_val, y_train, y_val = train_test_split(x_train_, y_train_, 
            #                                                   test_size=1/9, random_state = random_seed, stratify = y_train_)

            smote = SMOTE(random_state=random_seed)
            x_train, y_train = smote.fit_resample(x_train,y_train)
            y_train = y_train.reshape(-1,1)                              

            train_dataset = CustomDataset(x_train,y_train)
            val_dataset = CustomDataset(x_val,y_val)
            test_dataset = CustomDataset(x_test,y_test) 

            train_loader = DataLoader(dataset = train_dataset, batch_size = 64, shuffle = True)
            val_loader = DataLoader(dataset = val_dataset, batch_size = 64, shuffle = False)
            test_loader = DataLoader(dataset = test_dataset, batch_size = 64, shuffle = False)            


            for  lr in lr_list:
                for num_fc in num_fc_list:
                    val_auc_list = []
                    test_auc_list = []
                    self.model = PINNet(input_dim,pathway_info,num_fc)
                    self.optimizer = torch.optim.Adam(self.model.parameters(), lr=lr, weight_decay = 0)
                    self.criterion = nn.CrossEntropyLoss()
                    self.model = self.model.to(self.device)
                    early_stopping = EarlyStopping(patience=20, verbose = True, path = 'checkpoint_ES.pt')

                    ##train 
                    for epoch in range(0, 200):
                        #print(epoch)
                        for batch_idx, samples in enumerate(train_loader):
                            _,_ = self.train_step(samples,training = True)
                        ##early stopping
                        y_prob, y_true = [],[]
                        for batch_idx, samples in enumerate(val_loader):
                            prob, true = self.train_step(samples,training = False)

                            y_prob.extend(prob.detach().cpu().numpy())
                            y_true.extend(true.cpu().numpy())

                        val_auc, _, _, _ = self.evalutaion(y_true,y_prob)
                        #print(val_auc)

                        early_stopping(val_auc, self.model, epoch)
                        if early_stopping.early_stop:
                            break

                    ##validation 
                    self.model = torch.load('checkpoint_ES.pt')
                    y_prob, y_true = [],[]
                    for batch_idx, samples in enumerate(val_loader):     
                        prob, true = self.train_step(samples,training = False)
                        y_prob.extend(prob.detach().cpu().numpy())
                        y_true.extend(true.cpu().numpy())

                    val_auc, val_precision, val_recall, val_f1 = self.evalutaion(y_true,y_prob)

                    ##test
                    y_prob, y_true = [],[]
                    for batch_idx, samples in enumerate(test_loader):     
                        prob, true = self.train_step(samples,training = False)

                        y_prob.extend(prob.detach().cpu().numpy())
                        y_true.extend(true.cpu().numpy())
                    test_auc,_,_,_ = self.evalutaion(y_true,y_prob)
                    if test_auc > best_test_auc:
                        best_test_auc = test_auc
                        best_fold = fold+1
                        print("test auc of this fold: ", test_auc)
                        print("Best fold so far: ",best_fold)
                        ## train
                        y_prob_train, y_true_train = [],[]
                        for batch_idx, samples in enumerate(train_loader):     
                            prob, true = self.train_step(samples,training = False)
                            y_prob_train.extend(prob.detach().cpu().numpy())
                            y_true_train.extend(true.cpu().numpy())
                        ## independent test
                        y_prob_ind, y_true_ind = [], []
                        for batch_idx, samples in enumerate(test_loader1):     
                            prob, true = self.train_step(samples,training = False)
                            y_prob_ind.extend(prob.detach().cpu().numpy())
                            y_true_ind.extend(true.cpu().numpy())
                        shap_values, shap_values_feat = self.get_shap_values(x_train, expression1)
                    else:
                        print("test auc of this fold: ", test_auc)
                        print("Best fold so far: ",best_fold)
 
        return y_prob_ind, y_true_ind, y_prob_train, y_true_train, shap_values, shap_values_feat
        
    def train_step(self, batch_item, training):
        data,label = batch_item
        data = data.to(self.device)
        label = label.to(self.device)
        if training is True:
            self.model.train()
            self.optimizer.zero_grad()
            with torch.cuda.amp.autocast():
                out = self.model(data)
                true = torch.reshape(label,(-1,))
                loss = self.criterion(out,true)
                prob = out[:,1]                
                
            loss.backward()
            self.optimizer.step()
            return prob, true
        else:
            self.model.eval()
            with torch.no_grad():
                out = self.model(data)
                true = torch.reshape(label,(-1,))
                prob = out[:,1]   
            return prob, true  
        
    def evalutaion(self, y_true, y_prob):
        np.seterr(divide='ignore', invalid='ignore')
        auc = roc_auc_score(y_true,y_prob)
        precision,recall,_ = precision_recall_curve(y_true,y_prob)
        f1 = (2*precision*recall)/(precision+recall)
        idx = np.nanargmax(f1)
        pr = precision[idx] 
        rc = recall[idx] 
        f1 = f1[idx] 
        return auc, pr, rc, f1

    def seed_worker(self, random_seed):
        torch.manual_seed(random_seed)
        torch.cuda.manual_seed(random_seed)
        torch.cuda.manual_seed_all(random_seed)  # if use multi-GPU
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
        np.random.seed(random_seed)
        random.seed(random_seed)
        
    def get_shap_values(self, x_train, x_val):
            self.model.eval()
            explainer = shap.GradientExplainer((self.model, self.model.fc2), torch.tensor(x_train, dtype=torch.float32).to(self.device))
            shap_values = explainer.shap_values(torch.tensor(x_val, dtype=torch.float32).to(self.device))
            explainer_feat = shap.GradientExplainer(self.model, torch.tensor(x_train, dtype=torch.float32).to(self.device))
            shap_values_feat = explainer_feat.shap_values(torch.tensor(x_val, dtype=torch.float32).to(self.device))
            return shap_values, shap_values_feat
        
def independent_test():
    parser = argparse.ArgumentParser(description='Independent Test')
    parser.add_argument('--input_train_csv', default='gene.average.beta.by.intensity.train.csv', help='Input training csv file of X_i. Rows are genes, columns are samples.')
    parser.add_argument('--input_test_csv', default='gene.average.beta.by.intensity.test.csv', help='Input testing csv file of X_i. Rows are genes, columns are samples.')
    parser.add_argument('--input_label_train', default='label_train.csv', help='Input training label file path.')
    parser.add_argument('--input_label_test', default='label_test.csv', help='Input testing label file path.')
    parser.add_argument('--input_pathway', default='pathway_gobp.csv', help='Pathway information from dataprep.')
    parser.add_argument('--output', default='pred.pkl', help='Prediction output file path. pkl file.')
    parser.add_argument('--output_shap', default='shap.pkl', help='Model interpretation output file path. pkl file.')
    parser.add_argument('--mlp', action='store_true', help='DNN only model.')

    args = parser.parse_args()

    input_file_train_csv = args.input_train_csv
    input_file_test_csv = args.input_test_csv
    input_file_label_train = args.input_label_train
    input_file_label_test = args.input_label_test
    input_file_path = args.input_pathway
    output_file = args.output
    output_file_shap = args.output_shap
    mlp = args.mlp
    
    pathway_info, expression, status, features = preprocessing1(input_file_path, input_file_train_csv, input_file_label_train)
    pathway_info1, expression1, status1, features1 = preprocessing1(input_file_path, input_file_test_csv, input_file_label_test)
    
    trainArgs = {}
    trainArgs['pathway_info'] = pathway_info
    trainArgs['features'] = features
    trainArgs['num_fc_list'] = [64]
    trainArgs['lr_list'] = [0.0005]
    #trainArgs['num_fc_list'] = [32]
    #trainArgs['lr_list'] = [0.0001]
    trainArgs['device'] = '0'
    trainArgs['seed'] = 0
    trainArgs['pathway'] = pathway
    trainArgs['tissue'] = tissue
    trainArgs['filename'] = 'result.csv'
    
    train = train_test(trainArgs)
    y_prob, y_test, y_prob_train, y_true_train, shap_values, shap_values_feat= train.kfold()
    
    scores = []
    index = []
    
    y_true = np.array(y_test).astype(int)
    y_prob = np.array(y_prob)
    auc = roc_auc_score(y_true,y_prob)
    prauc = average_precision_score(y_true,y_prob)
    pr1, rc1, f1 = fscore(y_true,y_prob,1)
    pr2, rc2, f2 = fscore(y_true,y_prob,2)
    pr5, rc5, f5 = fscore(y_true,y_prob,5)
    pr10, rc10, f10 = fscore(y_true,y_prob,10)
    threshold = get_threshold(y_true_train_mlp, y_prob_train_mlp)
    mcc, threshold = mcc_score(y_true, y_prob,threshold)
    y_pred = np.greater(y_prob, threshold).astype(int)
    acc = accuracy_score(y_true, y_pred)
    tpr, tnr, ppv, npv = calculate_metrics(y_true, y_pred)
    scores.append([auc,prauc, mcc,threshold, acc,tpr, tnr, ppv, npv,f1,f2,f5,f10])
    index.append("Test")
    y_true = np.array(y_true_train).astype(int)
    y_prob = np.array(y_prob_train)
    auc = roc_auc_score(y_true,y_prob)
    prauc = average_precision_score(y_true,y_prob)
    pr1, rc1, f1 = fscore(y_true,y_prob,1)
    pr2, rc2, f2 = fscore(y_true,y_prob,2)
    pr5, rc5, f5 = fscore(y_true,y_prob,5)
    pr10, rc10, f10 = fscore(y_true,y_prob,10)
    threshold = get_threshold(y_true_train_mlp, y_prob_train_mlp)
    mcc, threshold = mcc_score(y_true, y_prob,threshold)
    y_pred = np.greater(y_prob, threshold).astype(int)
    acc = accuracy_score(y_true, y_pred)
    tpr, tnr, ppv, npv = calculate_metrics(y_true, y_pred)
    scores.append([auc,prauc, mcc,threshold, acc,tpr, tnr, ppv, npv,f1,f2,f5,f10])
    index.append("Train")

    shap = shap_values
    feature_importance = np.sum(np.abs(shap[1]), axis=0)[64:]
    print(feature_importance.shape)
    positive_sum = np.sum(shap[1] * (shap[1] > 0), axis=0)[64:]
    negative_sum = np.sum(shap[1] * (shap[1] < 0), axis=0)[64:]
    fi_result = pd.DataFrame({"path_name":features,"importance":feature_importance, "pos":positive_sum,"neg": np.abs(negative_sum)})
    df_sorted = fi_result.sort_values(by='importance', ascending=False)
    df_sorted.to_csv("Independent_Importance.csv",index=False)
    features = pd.read_csv("gene.average.beta.by.intensity.csv",index_col=0).index.tolist()
    print(len(features))
    shap = shaps[1]
    feature_importance = np.sum(np.abs(shap[1]), axis=0)
    print(feature_importance.shape)
    positive_sum = np.sum(shap[1] * (shap[1] > 0), axis=0)
    negative_sum = np.sum(shap[1] * (shap[1] < 0), axis=0)
    fi_result = pd.DataFrame({"path_name":features,"importance":feature_importance, "pos":positive_sum,"neg": np.abs(negative_sum)})
    df_sorted = fi_result.sort_values(by='importance', ascending=False)
    df_sorted.to_csv("Independent_FEP_GOBP_Feature_Importance.csv",index=False)