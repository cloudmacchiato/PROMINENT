import os, copy, torch, random, time, datetime
import numpy as np
import pandas as pd
from torch import nn
import torch.nn.functional as F
from torch.utils.data import Dataset
from torch.utils.data import DataLoader
from sklearn import metrics
from sklearn.feature_selection import mutual_info_classif
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.metrics import precision_recall_curve
from statistics import mean 
from imblearn.over_sampling import SMOTE
import shap
import pickle
from PROMINENT_methylation.model import *
from PROMINENT_methylation.utils import *
import matplotlib.pyplot as plt

class train_kfold:
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
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        directory = './CV_best_model'
        if not os.path.exists(directory):
            os.makedirs(directory)

    def kfold(self):
        trainArgs = self.trainArgs
        x_data = trainArgs['x_data']
        input_dim = x_data.shape[1]
        y_data = trainArgs['y_data']

        pathway_info = trainArgs['pathway_info'].to(self.device)
        num_pathway = pathway_info.shape[0]
        num_fc_list = trainArgs['num_fc_list']
        lr_list = trainArgs['lr_list']
        random_seed = trainArgs['seed']

        kfold = StratifiedKFold(n_splits = 5, shuffle=True, random_state = random_seed)
        
        result = pd.DataFrame(columns=['hyperparam','Fold', 'Valid_AUC','Valid_Precision','Valid_Recall','Valid_F1',
                                       'Test_AUC','Test_Precision','Test_Recall','Test_F1'])   

        for fold, (train_index, test_index) in enumerate(kfold.split(x_data, y_data)):   
            print('****************************************************************************')
            print('Fold {} / {}'.format(fold + 1 , kfold.get_n_splits()))
            print('****************************************************************************')
            x_train_ = x_data[train_index]
            y_train_ = y_data[train_index] 
            x_test = x_data[test_index]  
            y_test = y_data[test_index] 
            
            dict_data = {"train": train_index, "test": test_index}
            with open(f"../../data/folds/fold_{fold+1}.pkl", 'wb') as file:
                pickle.dump(dict_data, file)
            test_dataset = CustomDataset(x_test,y_test)    
            test_loader = DataLoader(dataset = test_dataset, batch_size = 64, shuffle = False)                   
            best_val_auc = 0
            for  lr in lr_list:
                for num_fc in num_fc_list:
                    val_auc_list = []
                    val_precision_list = []
                    val_recall_list = []
                    val_f1_list = []
                    seeds = [1,2,3,4,5,6,7,8,9,10]
                    for train_val_seed in seeds:
                        x_train, x_val, y_train, y_val, index_train, index_val = train_test_split(x_train_, y_train_, train_index, test_size=1/9, random_state = train_val_seed, stratify = y_train_)
                        smote = SMOTE(random_state=random_seed)
                        x_train, y_train = smote.fit_resample(x_train,y_train)
                        y_train = y_train.reshape(-1,1) 
                        train_dataset = CustomDataset(x_train,y_train)
                        val_dataset = CustomDataset(x_val,y_val)
                        train_loader = DataLoader(dataset = train_dataset, batch_size = 64, shuffle = True)
                        val_loader = DataLoader(dataset = val_dataset, batch_size = 64, shuffle = False)
                        self.model = PINNet(input_dim,pathway_info,num_fc)
                        self.optimizer = torch.optim.Adam(self.model.parameters(), lr=lr, weight_decay = 0)
                        self.criterion = nn.CrossEntropyLoss()
                        self.model = self.model.to(self.device)
                        early_stopping = EarlyStopping(patience=20, verbose = True, path = 'checkpoint_ES.pt')

                        ##train 
                        for epoch in range(0, 200):
                            for batch_idx, samples in enumerate(train_loader):
                                _,_ = self.train_step(samples,training = True)
                            ##early stopping
                            y_prob, y_true = [],[]
                            for batch_idx, samples in enumerate(val_loader):
                                prob, true = self.train_step(samples,training = False)

                                y_prob.extend(prob.detach().cpu().numpy())
                                y_true.extend(true.cpu().numpy())

                            val_auc, _, _, _ = self.evalutaion(y_true,y_prob)
                            print(val_auc)
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
                        val_auc_list.append(val_auc)
                        val_precision_list.append(val_precision)
                        val_recall_list.append(val_recall)
                        val_f1_list.append(val_f1)
                    mean_val_auc = np.mean(val_auc_list)
                    mean_val_precision = np.mean(val_precision_list)
                    mean_val_recall = np.mean(val_recall_list)
                    mean_val_f1 = np.mean(val_f1_list)
                    
                    if val_auc > best_val_auc:
                        best_val_auc = val_auc 

                    ##test
                    y_prob, y_true = [],[]
                    for batch_idx, samples in enumerate(test_loader):     
                        prob, true = self.train_step(samples,training = False)

                        y_prob.extend(prob.detach().cpu().numpy())
                        y_true.extend(true.cpu().numpy())

                    test_auc, test_precision, test_recall, test_f1 = self.evalutaion(y_true,y_prob)
                    
                    result = pd.concat([result, pd.DataFrame({'hyperparam': ["lr:{} / num_fc:{}".format(str(lr),str(num_fc))],'Fold':[fold],
                                            'Valid_AUC': [mean_val_auc], 'Valid_Precision': [mean_val_precision], 
                                            'Valid_Recall': [mean_val_recall], 'Valid_F1': [mean_val_f1],
                                            'Test_AUC': [test_auc], 'Test_Precision': [test_precision], 
                                            'Test_Recall': [test_recall], 'Test_F1': [test_f1]})], ignore_index=True)
        
        return result
        
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
        

class train_kfold1:
    def __init__(self, trainArgs):
   
        """
        args:
            x_data          : gene expression data
            y_data          : patient label
            pathway_info    : pathway matrics
            features        : feature names list
            num_fc_list     : number of fully connected nodes 
            lr_list         : learning rate
            device          : GPU device
        Returns:
            AUC, Precision, Recall, F1
        """
        self.trainArgs = trainArgs
        self.seed_worker(trainArgs['seed'])
        os.environ["CUDA_VISIBLE_DEVICES"]=trainArgs['device']
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        directory = './CV_best_model'
        if not os.path.exists(directory):
            os.makedirs(directory)

    def kfold(self):
        trainArgs = self.trainArgs
        x_data = trainArgs['x_data']
        input_dim = x_data.shape[1]
        y_data = trainArgs['y_data']

        pathway_info = trainArgs['pathway_info'].to(self.device)
        num_pathway = pathway_info.shape[0]
        num_fc_list = trainArgs['num_fc_list']
        lr_list = trainArgs['lr_list']
        random_seed = trainArgs['seed']

        kfold = StratifiedKFold(n_splits = 5, shuffle=True, random_state = random_seed)
        
        result = pd.DataFrame(columns=['hyperparam','Fold', 'Valid_AUC','Valid_Precision','Valid_Recall','Valid_F1',
                                       'Test_AUC','Test_Precision','Test_Recall','Test_F1'])

        for fold, (train_index, test_index) in enumerate(kfold.split(x_data, y_data)):   
            print('****************************************************************************')
            print('Fold {} / {}'.format(fold + 1 , kfold.get_n_splits()))
            print('****************************************************************************')
            x_train_ = x_data[train_index]
            y_train_ = y_data[train_index] 
            x_test = x_data[test_index]  
            y_test = y_data[test_index] 
            x_train, x_val, y_train, y_val = train_test_split(x_train_, y_train_, 
                                                              test_size=1/9, random_state = random_seed, stratify = y_train_)
            
            smote = SMOTE(random_state=random_seed)
            x_train, y_train = smote.fit_resample(x_train,y_train)
            y_train = y_train.reshape(-1,1)                              
            print(torch.tensor(x_train).shape)
            
            train_dataset = CustomDataset(x_train,y_train)
            val_dataset = CustomDataset(x_val,y_val)
            test_dataset = CustomDataset(x_test,y_test) 
                       
            train_loader = DataLoader(dataset = train_dataset, batch_size = 64, shuffle = True)
            val_loader = DataLoader(dataset = val_dataset, batch_size = 64, shuffle = False)
            test_loader = DataLoader(dataset = test_dataset, batch_size = 64, shuffle = False)            
                   
            best_val_auc = 0
            for  lr in lr_list:
                for num_fc in num_fc_list:
                    val_auc_list = []
                    test_auc_list = []
                    self.model = PINNet(input_dim,pathway_info,num_fc)
                    self.optimizer = torch.optim.Adam(self.model.parameters(), lr=lr, weight_decay = 0)
                    self.criterion = nn.CrossEntropyLoss()
                    self.model = self.model.to(self.device)
                    early_stopping = EarlyStopping(patience=10, verbose = True, path = 'checkpoint_ES.pt')
                   
                    ##train 
                    for epoch in range(0, 200):
                        for batch_idx, samples in enumerate(train_loader):
                            _,_ = self.train_step(samples,training = True)
                        ##early stopping
                        y_prob, y_true = [],[]
                        for batch_idx, samples in enumerate(val_loader):
                            prob, true = self.train_step(samples,training = False)

                            y_prob.extend(prob.detach().cpu().numpy())
                            y_true.extend(true.cpu().numpy())
                        val_auc, _, _, _ = self.evalutaion(y_true,y_prob)

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
                    print(">>[val] auc : {:.4f}, precision : {:.4f}, recall : {:.4f}, f1 : {:.4f}"
                          .format(val_auc, val_precision, val_recall, val_f1))                    
                    
                    if val_auc > best_val_auc:
                        best_val_auc = val_auc 
                    torch.save({'model': self.model.state_dict(), 'fold': fold,'learning_rate':lr,'num_fc':num_fc},
                               './CV_best_model/best_model_fold_{}.pt'.format(fold))
                    
                    ## Get SHAP values
                    features = trainArgs['features']
                    shap_values = self.get_shap_values(x_train, x_val)
                    feature_importance = pd.DataFrame({"Feature": features, f'Importance{fold+1}': np.mean(np.abs(shap_values), axis=0)})
                    top20 = feature_importance.sort_values(by=f'Importance{fold+1}').head(20)
                    plt.barh(top20['Feature'], top20[f'Importance{fold+1}'])
                    plt.title(f'Top 20 Features by Importance of Fold {fold+1}')
                    plt.xlabel('Importance')
                    plt.ylabel('Feature')
                    plt.show()
                    if fold == 0:
                        fi_result = feature_importance
                    else:
                        fi_result = pd.concat([fi_result, feature_importance.iloc[:,1]],axis=1)
                    
                    ##test
                    y_prob, y_true = [],[]
                    for batch_idx, samples in enumerate(test_loader):     
                        prob, true = self.train_step(samples,training = False)

                        y_prob.extend(prob.detach().cpu().numpy())
                        y_true.extend(true.cpu().numpy())

                    test_auc, test_precision, test_recall, test_f1 = self.evalutaion(y_true,y_prob)

                    print(">>[test] auc : {:.4f}, precision : {:.4f}, recall : {:.4f}, f1 : {:.4f}\n"
                          .format(test_auc, test_precision, test_recall, test_f1))
                    
                    result = pd.concat([result, pd.DataFrame({'hyperparam': ["lr:{} / num_fc:{}".format(str(lr),str(num_fc))],'Fold':[fold],
                                            'Valid_AUC': [val_auc], 'Valid_Precision': [val_precision], 
                                            'Valid_Recall': [val_recall], 'Valid_F1': [val_f1],
                                            'Test_AUC': [test_auc], 'Test_Precision': [test_precision], 
                                            'Test_Recall': [test_recall], 'Test_F1': [test_f1]})], ignore_index=True)
        
        return result, fi_result
        
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
            e = shap.DeepExplainer(self.model, data)
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
        explainer = shap.GradientExplainer(self.model, torch.tensor(x_train, dtype=torch.float32).to(self.device))
        shap_values = explainer.shap_values(torch.tensor(x_val, dtype=torch.float32).to(self.device))[1]
        return np.array(shap_values)
    
class train_kfold2:
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
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        print(self.device)

    def kfold(self):
        print(self.device)
        trainArgs = self.trainArgs
        x_data = trainArgs['x_data']
        input_dim = x_data.shape[1]
        y_data = trainArgs['y_data']

        pathway_info = trainArgs['pathway_info'].to(self.device)
        num_pathway = pathway_info.shape[0]
        num_fc_list = trainArgs['num_fc_list']
        lr_list = trainArgs['lr_list']
        random_seed = trainArgs['seed']

        kfold = StratifiedKFold(n_splits = 5, shuffle=True, random_state = random_seed)
        
        result = pd.DataFrame(columns=['hyperparam','Fold', 'Valid_AUC','Valid_Precision','Valid_Recall','Valid_F1',
                                       'Test_AUC','Test_Precision','Test_Recall','Test_F1'])   
        shap_ls = []
        shap_feat_ls = []
        x_tests = []
        test_preds = []
        preds_fn = trainArgs['filename']
        shap_fn = trainArgs['filename2']
        for fold, (train_index, test_index) in enumerate(kfold.split(x_data, y_data)):   
            print('****************************************************************************')
            print('Fold {} / {}'.format(fold + 1 , kfold.get_n_splits()))
            print('****************************************************************************')
            x_train_ = x_data[train_index]
            y_train_ = y_data[train_index] 
            x_test = x_data[test_index]  
            y_test = y_data[test_index] 
            x_train, x_val, y_train, y_val = train_test_split(x_train_, y_train_, 
                                                              test_size=1/9, random_state = random_seed, stratify = y_train_)
            dict_data = {"train": train_index, "test": test_index}
            
            smote = SMOTE(random_state=random_seed)
            x_train, y_train = smote.fit_resample(x_train,y_train)
            y_train = y_train.reshape(-1,1)  
            
            train_dataset = CustomDataset(x_train,y_train)
            val_dataset = CustomDataset(x_val,y_val)
            test_dataset = CustomDataset(x_test,y_test)  
                       
            train_loader = DataLoader(dataset = train_dataset, batch_size = 64, shuffle = True)
            val_loader = DataLoader(dataset = val_dataset, batch_size = 64, shuffle = False)
            test_loader = DataLoader(dataset = test_dataset, batch_size = 64, shuffle = False)            
                   
            best_val_auc = 0
            for  lr in lr_list:
                for num_fc in num_fc_list:
                    val_auc_list = []
                    test_auc_list = []
                    self.model = PINNet(input_dim,pathway_info,num_fc)
                    self.optimizer = torch.optim.Adam(self.model.parameters(), lr=lr, weight_decay = 0)
                    self.criterion = nn.CrossEntropyLoss()
                    self.model = self.model.to(self.device)
                    early_stopping = EarlyStopping(patience=10, verbose = True, path = f'checkpoint_ES.pt')
                   
                    ##train 
                    for epoch in range(0, 200):
                        for batch_idx, samples in enumerate(train_loader):
                            _,_ = self.train_step(samples,training = True)
                        ##early stopping
                        y_prob, y_true = [],[]
                        for batch_idx, samples in enumerate(val_loader):
                            prob, true = self.train_step(samples,training = False)

                            y_prob.extend(prob.detach().cpu().numpy())
                            y_true.extend(true.cpu().numpy())

                        val_auc, _, _, _, _ = self.evalutaion(y_true,y_prob)

                        early_stopping(val_auc, self.model, epoch)
                        if early_stopping.early_stop:
                            break
                            
                    ##validation 
                    self.model = torch.load(f'checkpoint_ES.pt')
                    y_prob, y_true = [],[]
                    for batch_idx, samples in enumerate(val_loader):     
                        prob, true = self.train_step(samples,training = False)
                        y_prob.extend(prob.detach().cpu().numpy())
                        y_true.extend(true.cpu().numpy())

                    val_auc, val_precision, val_recall, val_f1, _ = self.evalutaion(y_true,y_prob)
                 
                    
                    if val_auc > best_val_auc:
                        best_val_auc = val_auc 

                    ##test
                    y_prob, y_true = [],[]
                    for batch_idx, samples in enumerate(test_loader):     
                        prob, true = self.train_step(samples,training = False)

                        y_prob.extend(prob.detach().cpu().numpy())
                        y_true.extend(true.cpu().numpy())

                    test_auc, test_precision, test_recall, test_f1, test_pr_auc = self.evalutaion(y_true,y_prob)
                    test_preds.append([y_true,y_prob])
                    result = pd.concat([result, pd.DataFrame({'hyperparam': ["lr:{} / num_fc:{}".format(str(lr),str(num_fc))],'Fold':[fold],
                                            'Valid_AUC': [val_auc], 'Valid_Precision': [val_precision], 
                                            'Valid_Recall': [val_recall], 'Valid_F1': [val_f1],
                                            'Test_AUC': [test_auc], 'Test_Precision': [test_precision], 
                                            'Test_Recall': [test_recall], 'Test_F1': [test_f1], 'Test_PrAUC': [test_pr_auc]})], ignore_index=True)
                    
                    ##SHAP
                    shap_values, shap_values_feat = self.get_shap_values(x_train, x_test)
                    shap_ls.append(shap_values)
                    shap_feat_ls.append(shap_values_feat)
                    x_tests.append(x_test)
        with open(shap_fn, 'wb') as file:
            pickle.dump({"path":shap_ls,"feat":shap_feat_ls}, file)
        with open(preds_fn, 'wb') as file:
            pickle.dump(test_preds, file)
        return result
        
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
        pr_auc = average_precision_score(y_true,y_prob)
        precision,recall,_ = precision_recall_curve(y_true,y_prob)
        f1 = (2*precision*recall)/(precision+recall)
        idx = np.nanargmax(f1)
        pr = precision[idx] 
        rc = recall[idx] 
        f1 = f1[idx] 
        return auc, pr, rc, f1, pr_auc

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

class train_kfold_mlp:
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
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        directory = './CV_best_model'
        if not os.path.exists(directory):
            os.makedirs(directory)

    def kfold(self):
        trainArgs = self.trainArgs
        x_data = trainArgs['x_data']
        input_dim = x_data.shape[1]
        y_data = trainArgs['y_data']
        pathway_info = trainArgs['pathway_info'].to(self.device)
        num_pathway = pathway_info.shape[0]
        num_fc_list = trainArgs['num_fc_list']
        lr_list = trainArgs['lr_list']
        random_seed = trainArgs['seed']

        kfold = StratifiedKFold(n_splits = 5, shuffle=True, random_state = random_seed)
        
        result = pd.DataFrame(columns=['hyperparam','Fold', 'Valid_AUC','Valid_Precision','Valid_Recall','Valid_F1',
                                       'Test_AUC','Test_Precision','Test_Recall','Test_F1'])   
        
        test_preds = []
        preds_fn = trainArgs['filename']
        for fold, (train_index, test_index) in enumerate(kfold.split(x_data, y_data)):   
            print('****************************************************************************')
            print('Fold {} / {}'.format(fold + 1 , kfold.get_n_splits()))
            print('****************************************************************************')
            x_train_ = x_data[train_index]
            y_train_ = y_data[train_index] 
            x_test = x_data[test_index]  
            y_test = y_data[test_index] 
            x_train, x_val, y_train, y_val = train_test_split(x_train_, y_train_, 
                                                              test_size=1/9, random_state = random_seed, stratify = y_train_)
            
            smote = SMOTE(random_state=random_seed)
            x_train, y_train = smote.fit_resample(x_train,y_train)
            y_train = y_train.reshape(-1,1)                              
            
 
            train_dataset = CustomDataset(x_train,y_train)
            val_dataset = CustomDataset(x_val,y_val)
            test_dataset = CustomDataset(x_test,y_test)
                       
            train_loader = DataLoader(dataset = train_dataset, batch_size = 64, shuffle = True)
            val_loader = DataLoader(dataset = val_dataset, batch_size = 64, shuffle = False)
            test_loader = DataLoader(dataset = test_dataset, batch_size = 64, shuffle = False)            
                   
            best_val_auc = 0
            for  lr in lr_list:
                for num_fc in num_fc_list:
                    val_auc_list = []
                    test_auc_list = []
                    self.model = MLP(input_dim,pathway_info,num_fc)
                    self.optimizer = torch.optim.Adam(self.model.parameters(), lr=lr, weight_decay = 0)
                    self.criterion = nn.CrossEntropyLoss()
                    self.model = self.model.to(self.device)
                    early_stopping = EarlyStopping(patience=20, verbose = True, path = f'checkpoint_ES.pt')
                   
                    ##train 
                    for epoch in range(0, 200):
                        for batch_idx, samples in enumerate(train_loader):
                            _,_ = self.train_step(samples,training = True)
                        ##early stopping
                        y_prob, y_true = [],[]
                        for batch_idx, samples in enumerate(val_loader):
                            prob, true = self.train_step(samples,training = False)

                            y_prob.extend(prob.detach().cpu().numpy())
                            y_true.extend(true.cpu().numpy())

                        val_auc, _, _, _, _ = self.evalutaion(y_true,y_prob)

                        early_stopping(val_auc, self.model, epoch)
                        if early_stopping.early_stop:
                            break
                            
                    ##validation 
                    self.model = torch.load(f'checkpoint_ES.pt')
                    y_prob, y_true = [],[]
                    for batch_idx, samples in enumerate(val_loader):     
                        prob, true = self.train_step(samples,training = False)
                        y_prob.extend(prob.detach().cpu().numpy())
                        y_true.extend(true.cpu().numpy())

                    val_auc, val_precision, val_recall, val_f1, _ = self.evalutaion(y_true,y_prob)
                 
                    
                    if val_auc > best_val_auc:
                        best_val_auc = val_auc 

                    ##test
                    y_prob, y_true = [],[]
                    for batch_idx, samples in enumerate(test_loader):     
                        prob, true = self.train_step(samples,training = False)

                        y_prob.extend(prob.detach().cpu().numpy())
                        y_true.extend(true.cpu().numpy())

                    test_auc, test_precision, test_recall, test_f1, test_pr_auc = self.evalutaion(y_true,y_prob)
                    test_preds.append([y_true,y_prob])
                    result = pd.concat([result, pd.DataFrame({'hyperparam': ["lr:{} / num_fc:{}".format(str(lr),str(num_fc))],'Fold':[fold],
                                            'Valid_AUC': [val_auc], 'Valid_Precision': [val_precision], 
                                            'Valid_Recall': [val_recall], 'Valid_F1': [val_f1],
                                            'Test_AUC': [test_auc], 'Test_Precision': [test_precision], 
                                            'Test_Recall': [test_recall], 'Test_F1': [test_f1], 'Test_PrAUC': [test_pr_auc]})], ignore_index=True)
        with open(preds_fn, 'wb') as file:
            pickle.dump(test_preds, file)
        return result
    
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
        pr_auc = average_precision_score(y_true,y_prob)
        precision,recall,_ = precision_recall_curve(y_true,y_prob)
        f1 = (2*precision*recall)/(precision+recall)
        idx = np.nanargmax(f1)
        pr = precision[idx] 
        rc = recall[idx] 
        f1 = f1[idx] 
        return auc, pr, rc, f1, pr_auc

    def seed_worker(self, random_seed):
        torch.manual_seed(random_seed)
        torch.cuda.manual_seed(random_seed)
        torch.cuda.manual_seed_all(random_seed)  # if use multi-GPU
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
        np.random.seed(random_seed)
        random.seed(random_seed)
               
class train_kfold3:
    def __init__(self, trainArgs):
   
        """
        args:
            x_data          : Matrix data with gene features and methylation
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
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    def kfold(self):
        trainArgs = self.trainArgs
        x_data = trainArgs['x_data']
        y_data = trainArgs['y_data']
        pathway_info = trainArgs['pathway_info']
        sel_feat_num = trainArgs['sel_feat_num']
        num_fc_list = trainArgs['num_fc_list']
        lr_list = trainArgs['lr_list']
        random_seed = trainArgs['seed']

        kfold = StratifiedKFold(n_splits = 5, shuffle=True, random_state = random_seed)
        
        result = pd.DataFrame(columns=['hyperparam','Fold', 'Valid_AUC','Valid_Precision','Valid_Recall','Valid_F1',
                                       'Test_AUC','Test_Precision','Test_Recall','Test_F1'])   
        shap_ls = []
        x_tests = []
        sel_feat_idx = {}
        test_preds = []
        preds_fn = trainArgs['filename']
        feat_idx_fn = trainArgs['sel_feat_fn']
        for fold, (train_index, test_index) in enumerate(kfold.split(x_data, y_data)):   
            print('****************************************************************************')
            print('Fold {} / {}'.format(fold + 1 , kfold.get_n_splits()))
            print('****************************************************************************')
            x_train_ = x_data[train_index]
            y_train_ = y_data[train_index] 
            x_test = x_data[test_index]  
            y_test = y_data[test_index] 
            x_train, x_val, y_train, y_val = train_test_split(x_train_, y_train_, 
                                                              test_size=1/9, random_state = random_seed, stratify = y_train_)
            
            smote = SMOTE(random_state=random_seed)
            x_train, y_train = smote.fit_resample(x_train,y_train)
            y_train = y_train.reshape(-1,1)                              

            if os.path.exists(feat_idx_fn):
                with open(feat_idx_fn, 'rb') as file:
                    sel_feat_idx = pickle.load(file)
                    kept_indice1 = sel_feat_idx[f'fold{fold}'][0]
                    kept_indice2 = sel_feat_idx[f'fold{fold}'][1]
            else:
                print("Run feature selection first!")
            
            pathway_info_input = pathway_info[:,kept_indice2]
            pathway_info_input = pathway_info_input.to(self.device)
            input_size_gene = len(kept_indice2)
            input_size_meth = sel_feat_num
            
            train_dataset = CustomDataset(x_train[:,kept_indice1],y_train)
            val_dataset = CustomDataset(x_val[:,kept_indice1],y_val)
            test_dataset = CustomDataset(x_test[:,kept_indice1],y_test) 
                       
            train_loader = DataLoader(dataset = train_dataset, batch_size = 64, shuffle = True)
            val_loader = DataLoader(dataset = val_dataset, batch_size = 64, shuffle = False)
            test_loader = DataLoader(dataset = test_dataset, batch_size = 64, shuffle = False)            
                   
            best_val_auc = 0
            for  lr in lr_list:
                for num_fc in num_fc_list:
                    val_auc_list = []
                    test_auc_list = []
                    self.model = PINNet3(input_size_gene,sel_feat_num,pathway_info_input,num_fc)
                    self.optimizer = torch.optim.Adam(self.model.parameters(), lr=lr, weight_decay = 0)
                    self.criterion = nn.CrossEntropyLoss()
                    self.model = self.model.to(self.device)
                    early_stopping = EarlyStopping(patience=10, verbose = True, path = f'/ix/hpark/laz64/ipf/dlmethylation/run{random_seed}/checkpoint_ES.pt')
                   
                    ##train 
                    for epoch in range(0, 200):
                        for batch_idx, samples in enumerate(train_loader):
                            #print(samples[0].shape)
                            _,_ = self.train_step(samples,training = True, num_gene = input_size_gene)
                        ##early stopping
                        y_prob, y_true = [],[]
                        for batch_idx, samples in enumerate(val_loader):
                            prob, true = self.train_step(samples,training = False, num_gene = input_size_gene)
                            #print(true, prob)
                            y_prob.extend(prob.detach().cpu().numpy())
                            y_true.extend(true.cpu().numpy())
                        val_auc, _, _, _, _ = self.evalutaion(y_true,y_prob)

                        early_stopping(val_auc, self.model, epoch)
                        if early_stopping.early_stop:
                            break
                            
                    ##validation 
                    self.model = torch.load(f'/ix/hpark/laz64/ipf/dlmethylation/run{random_seed}/checkpoint_ES.pt')
                    y_prob, y_true = [],[]
                    for batch_idx, samples in enumerate(val_loader):     
                        prob, true = self.train_step(samples,training = False, num_gene = input_size_gene)
                        y_prob.extend(prob.detach().cpu().numpy())
                        y_true.extend(true.cpu().numpy())

                    val_auc, val_precision, val_recall, val_f1, _ = self.evalutaion(y_true,y_prob)
                 
                    
                    if val_auc > best_val_auc:
                        best_val_auc = val_auc 

                    ##test
                    y_prob, y_true = [],[]
                    for batch_idx, samples in enumerate(test_loader):     
                        prob, true = self.train_step(samples,training = False, num_gene = input_size_gene)

                        y_prob.extend(prob.detach().cpu().numpy())
                        y_true.extend(true.cpu().numpy())

                    test_auc, test_precision, test_recall, test_f1, test_pr_auc = self.evalutaion(y_true,y_prob)
                    test_preds.append([y_true,y_prob])
                    
                    result = pd.concat([result, pd.DataFrame({'hyperparam': ["lr:{} / num_fc:{}".format(str(lr),str(num_fc))],'Fold':[fold],
                                            'Valid_AUC': [val_auc], 'Valid_Precision': [val_precision], 
                                            'Valid_Recall': [val_recall], 'Valid_F1': [val_f1],
                                            'Test_AUC': [test_auc], 'Test_Precision': [test_precision], 
                                            'Test_Recall': [test_recall], 'Test_F1': [test_f1], 'Test_PrAUC': [test_pr_auc]})], ignore_index=True)
                    
                    ##SHAP
                    # shap_values = self.get_shap_values(x_train[:,kept_indice], x_test[:,kept_indice],num_gene=input_size_gene)
                    # shap_ls.append(shap_values)
                    # x_tests.append(x_test)
        # with open("shap.pkl", 'wb') as file:
        #     pickle.dump(shap_ls, file)
        with open(preds_fn, 'wb') as file:
            pickle.dump(test_preds, file)
        return result
        
    def train_step(self, batch_item, training,num_gene):
        data,label = batch_item
        data = data.to(self.device)
        input_gene = data[:,:num_gene]
        #print(input_gene.shape)
        input_meth = data[:,num_gene:]
        #print(input_meth.shape)
        label = label.to(self.device)
        if training is True:
            self.model.train()
            self.optimizer.zero_grad()
            with torch.cuda.amp.autocast():
                out = self.model(input_gene,input_meth)
                true = torch.reshape(label,(-1,))
                loss = self.criterion(out,true)
                prob = out[:,1]                
            loss.backward()
            self.optimizer.step()
            return prob, true
        else:
            self.model.eval()
            with torch.no_grad():
                print(input_gene)
                print(input_meth)
                out = self.model(input_gene,input_meth)
                print(out)
                true = torch.reshape(label,(-1,))
                prob = out[:,1]   
            return prob, true  
        
    def evalutaion(self, y_true, y_prob):
        np.seterr(divide='ignore', invalid='ignore')
        auc = roc_auc_score(y_true,y_prob)
        pr_auc = average_precision_score(y_true,y_prob)
        precision,recall,_ = precision_recall_curve(y_true,y_prob)
        f1 = (2*precision*recall)/(precision+recall)
        idx = np.nanargmax(f1)
        pr = precision[idx] 
        rc = recall[idx] 
        f1 = f1[idx] 
        return auc, pr, rc, f1, pr_auc

    def seed_worker(self, random_seed):
        torch.manual_seed(random_seed)
        torch.cuda.manual_seed(random_seed)
        torch.cuda.manual_seed_all(random_seed)  # if use multi-GPU
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
        np.random.seed(random_seed)
        random.seed(random_seed)
        
    def get_shap_values(self, x_train, x_val, num_gene):
            self.model.eval()
            input_gene = torch.tensor(x_train[:,:num_gene], dtype=torch.float32).to(self.device)
            input_meth = torch.tensor(x_train[:,num_gene:], dtype=torch.float32).to(self.device)
            input_gene_val = torch.tensor(x_val[:,:num_gene], dtype=torch.float32).to(self.device)
            input_meth_val = torch.tensor(x_val[:,num_gene:], dtype=torch.float32).to(self.device)
            explainer = shap.GradientExplainer((self.model, self.model.fc2), [input_gene,input_meth])
            shap_values = explainer.shap_values([input_gene_val, input_meth_val])
            return shap_values
        
class train_kfold4:
    def __init__(self, trainArgs):
   
        """
        args:
            x_data          : Matrix data with gene features and methylation
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
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        directory = './CV_best_model'
        if not os.path.exists(directory):
            os.makedirs(directory)

    def kfold(self):
        trainArgs = self.trainArgs
        x_data = trainArgs['x_data']
        y_data = trainArgs['y_data']
        pathway_info = trainArgs['pathway_info']
        sel_feat_num = trainArgs['sel_feat_num']
        num_fc_list = trainArgs['num_fc_list']
        lr_list = trainArgs['lr_list']
        random_seed = trainArgs['seed']
        input_size_gene = pathway_info.shape[1]
        input_size_meth = sel_feat_num
        pathway_info_input = pathway_info.to(self.device)
        kfold = StratifiedKFold(n_splits = 5, shuffle=True, random_state = random_seed)
        
        result = pd.DataFrame(columns=['hyperparam','Fold', 'Valid_AUC','Valid_Precision','Valid_Recall','Valid_F1',
                                       'Test_AUC','Test_Precision','Test_Recall','Test_F1'])   
        shap_ls = []
        x_tests = []
        sel_feat_idx = {}
        test_preds = []
        for fold, (train_index, test_index) in enumerate(kfold.split(x_data, y_data)):   
            print('****************************************************************************')
            print('Fold {} / {}'.format(fold + 1 , kfold.get_n_splits()))
            print('****************************************************************************')
            x_train_ = x_data[train_index]
            y_train_ = y_data[train_index] 
            x_test = x_data[test_index]  
            y_test = y_data[test_index] 
            x_train, x_val, y_train, y_val = train_test_split(x_train_, y_train_, 
                                                              test_size=1/9, random_state = random_seed, stratify = y_train_)
            
            smote = SMOTE(random_state=random_seed)
            x_train, y_train = smote.fit_resample(x_train,y_train)
            y_train = y_train.reshape(-1,1)                              

            if os.path.exists("mi_feat_idx_gobp.pkl"):
                with open("mi_feat_idx_gobp.pkl", 'rb') as file:
                    sel_feat_idx = pickle.load(file)
                    kept_indice3 = sel_feat_idx[f'fold{fold}'][2]
            else:
                print("Run feature selection first!")
            
            kept_indice1 = [i for i in range(input_size_gene)] + [x + input_size_gene for x in kept_indice3]
            print(len(kept_indice1),input_size_gene)
            train_dataset = CustomDataset(x_train[:,kept_indice1],y_train)
            val_dataset = CustomDataset(x_val[:,kept_indice1],y_val)
            test_dataset = CustomDataset(x_test[:,kept_indice1],y_test) 
                       
            train_loader = DataLoader(dataset = train_dataset, batch_size = 64, shuffle = True)
            val_loader = DataLoader(dataset = val_dataset, batch_size = 64, shuffle = False)
            test_loader = DataLoader(dataset = test_dataset, batch_size = 64, shuffle = False)            
                   
            best_val_auc = 0
            for  lr in lr_list:
                for num_fc in num_fc_list:
                    val_auc_list = []
                    test_auc_list = []
                    self.model = PINNet3(input_size_gene,sel_feat_num,pathway_info_input,num_fc)
                    self.optimizer = torch.optim.Adam(self.model.parameters(), lr=lr, weight_decay = 0)
                    self.criterion = nn.CrossEntropyLoss()
                    self.model = self.model.to(self.device)
                    early_stopping = EarlyStopping(patience=10, verbose = True, path = 'checkpoint_ES.pt')
                   
                    ##train 
                    for epoch in range(0, 200):
                        for batch_idx, samples in enumerate(train_loader):
                            #print(samples[0].shape)
                            _,_ = self.train_step(samples,training = True, num_gene = input_size_gene)
                        ##early stopping
                        y_prob, y_true = [],[]
                        for batch_idx, samples in enumerate(val_loader):
                            #print(samples)
                            
                            prob, true = self.train_step(samples,training = False, num_gene = input_size_gene)
                            #print(true, prob)
                            y_prob.extend(prob.detach().cpu().numpy())
                            y_true.extend(true.cpu().numpy())
                        #print(y_true, y_prob)
                        val_auc, _, _, _, _ = self.evalutaion(y_true,y_prob)

                        early_stopping(val_auc, self.model, epoch)
                        if early_stopping.early_stop:
                            break
                            
                    ##validation 
                    self.model = torch.load('checkpoint_ES.pt')
                    y_prob, y_true = [],[]
                    for batch_idx, samples in enumerate(val_loader):     
                        prob, true = self.train_step(samples,training = False, num_gene = input_size_gene)
                        y_prob.extend(prob.detach().cpu().numpy())
                        y_true.extend(true.cpu().numpy())

                    val_auc, val_precision, val_recall, val_f1, _ = self.evalutaion(y_true,y_prob)
                 
                    
                    if val_auc > best_val_auc:
                        best_val_auc = val_auc 

                    ##test
                    y_prob, y_true = [],[]
                    for batch_idx, samples in enumerate(test_loader):     
                        prob, true = self.train_step(samples,training = False, num_gene = input_size_gene)

                        y_prob.extend(prob.detach().cpu().numpy())
                        y_true.extend(true.cpu().numpy())

                    test_auc, test_precision, test_recall, test_f1, test_pr_auc = self.evalutaion(y_true,y_prob)
                    
                    test_preds.append([y_true,y_prob])
                    result = pd.concat([result, pd.DataFrame({'hyperparam': ["lr:{} / num_fc:{}".format(str(lr),str(num_fc))],'Fold':[fold],
                                            'Valid_AUC': [val_auc], 'Valid_Precision': [val_precision], 
                                            'Valid_Recall': [val_recall], 'Valid_F1': [val_f1],
                                            'Test_AUC': [test_auc], 'Test_Precision': [test_precision], 
                                            'Test_Recall': [test_recall], 'Test_F1': [test_f1], 'Test_PrAUC': [test_pr_auc]})], ignore_index=True)
                    
                    ##SHAP
                    # shap_values = self.get_shap_values(x_train[:,kept_indice], x_test[:,kept_indice],num_gene=input_size_gene)
                    # shap_ls.append(shap_values)
                    # x_tests.append(x_test)
        # with open("shap.pkl", 'wb') as file:
        #     pickle.dump(shap_ls, file)
        with open("test_preds.pkl", 'wb') as file:
            pickle.dump(test_preds, file)
        return result
        
    def train_step(self, batch_item, training,num_gene):
        data,label = batch_item
        data = data.to(self.device)
        input_gene = data[:,:num_gene]
        #print(input_gene.shape)
        input_meth = data[:,num_gene:]
        #print(input_meth.shape)
        label = label.to(self.device)
        if training is True:
            self.model.train()
            self.optimizer.zero_grad()
            with torch.cuda.amp.autocast():
                out = self.model(input_gene,input_meth)
                true = torch.reshape(label,(-1,))
                loss = self.criterion(out,true)
                prob = out[:,1]                
            loss.backward()
            self.optimizer.step()
            return prob, true
        else:
            self.model.eval()
            with torch.no_grad():
                out = self.model(input_gene,input_meth)
                true = torch.reshape(label,(-1,))
                prob = out[:,1]   
            return prob, true  
        
    def evalutaion(self, y_true, y_prob):
        np.seterr(divide='ignore', invalid='ignore')
        auc = roc_auc_score(y_true,y_prob)
        pr_auc = average_precision_score(y_true,y_prob)
        precision,recall,_ = precision_recall_curve(y_true,y_prob)
        f1 = (2*precision*recall)/(precision+recall)
        idx = np.nanargmax(f1)
        pr = precision[idx] 
        rc = recall[idx] 
        f1 = f1[idx] 
        return auc, pr, rc, f1, pr_auc

    def seed_worker(self, random_seed):
        torch.manual_seed(random_seed)
        torch.cuda.manual_seed(random_seed)
        torch.cuda.manual_seed_all(random_seed)  # if use multi-GPU
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
        np.random.seed(random_seed)
        random.seed(random_seed)
        
    def get_shap_values(self, x_train, x_val, num_gene):
            self.model.eval()
            input_gene = torch.tensor(x_train[:,:num_gene], dtype=torch.float32).to(self.device)
            input_meth = torch.tensor(x_train[:,num_gene:], dtype=torch.float32).to(self.device)
            input_gene_val = torch.tensor(x_val[:,:num_gene], dtype=torch.float32).to(self.device)
            input_meth_val = torch.tensor(x_val[:,num_gene:], dtype=torch.float32).to(self.device)
            explainer = shap.GradientExplainer((self.model, self.model.fc2), [input_gene,input_meth])
            shap_values = explainer.shap_values([input_gene_val, input_meth_val])
            return shap_values