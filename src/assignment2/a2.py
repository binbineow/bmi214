#the code ran on cardinal.stanford.edu
path0 = '/afs/.ir/users/b/c/bchen45/biomedin214/bmi214/'
path_local = '/Users/binbineow2/Dropbox/Russ/python_code/Common_Work_Station/BMI214/'
file_name_data = 'leukemia.csv'
from sklearn.neighbors import KNeighborsClassifier
import numpy as np
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import confusion_matrix
from sklearn.naive_bayes import GaussianNB
from sklearn.metrics import classification_report


def get_data(path0,file_name_data):
    x_data = []
    y_label = []
    #0 = AML
    #1 = ALL
    for line0 in open(path0+file_name_data):
        if 'CCND3' in line0:
            line0 = line0.rstrip().split(',')
            gene_label = line0[:-1]
        elif len(line0)>10:
            line0 = line0.rstrip().split(',')
            x_data.append([int(w) for w in line0[0:-1]])
            if line0[-1] == 'AML':
                y_label.append(0)
            else:
                y_label.append(1)

    return [np.array(x_data),np.array(y_label),np.array(gene_label)]
    
def cal_parameters(c_mat):
    #Thus in binary classification, the count of true negatives
    # is C_{0,0}, false negatives is C_{1,0}, true positives is C_{1,1} and false positives is C_{0,1}.
    #0 is negative (AML) 1 is ALL
    #thus tp for all = c 1,1/ (c1,1 + c1,0)
    #fp for all c(0,1)/(c0,0 + c0,1)
    #tp for aml = 1 - fp for all
    #fp for aml = 1 - tp for all
    tp_all = c_mat[1,1]/float(c_mat[1,1]+c_mat[1,0])
    fp_all = c_mat[0,1]/float(c_mat[0,0]+c_mat[0,1])
    tp_aml = 1 - fp_all
    fp_aml = 1 - tp_all
    a0 = float(c_mat[1,1]+c_mat[0,0])/(c_mat[1,1]+c_mat[1,0]+c_mat[0,0]+c_mat[0,1])
    print('TP for AML')
    print(tp_aml)
    print('TP for ALL')
    print(tp_all)
    print('FP for AML')
    print(fp_aml)
    print('FP for ALL')
    print(fp_all)
    print('Accuracy')
    print(a0)
    
    
def loo_knn(x_data,y_label):
    from sklearn.neighbors import KNeighborsClassifier
    loo = LeaveOneOut()
    y_pred = []
    y_real = []
    for train_mask, test_mask in loo.split(x_data):
        nkk0 = KNeighborsClassifier(n_neighbors=40).fit(x_data[train_mask],y_label[train_mask])
        y_pred.append(nkk0.predict(x_data[test_mask]))
        y_real.append(y_label[test_mask])
    c_mat = confusion_matrix(y_real, y_pred)
    print("KNN")
    cal_parameters(c_mat)
    #print(len(y_real))
    #print(len(y_pred))
    #print(classification_report(list(y_real),list(y_pred)))
    
def loo_nb(x_data,y_label):
    loo = LeaveOneOut()
    y_pred = []
    y_real = []
    for train_mask, test_mask in loo.split(x_data):
        nb0 = GaussianNB()
        nb0.fit(x_data[train_mask],y_label[train_mask])
        y_pred.append(nb0.predict(x_data[test_mask]))
        y_real.append(y_label[test_mask])
    c_mat = confusion_matrix(y_real, y_pred)
    print("NB")
    cal_parameters(c_mat)  
    #classification_report(y_real,y_pred) 
        


def main():
    [x_data,y_label,gene_label] = get_data(path_local,file_name_data)
    print(len(x_data))
    loo_knn(x_data,y_label)
    loo_nb(x_data,y_label)
    
    

main()