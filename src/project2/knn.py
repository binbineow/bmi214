#!/usr/bin/env python
#BMI214 Project 2
#Binbin Chen
#bchen45

#For the KNN algorithm, I start with a training set of data points with class labels 
#ALL positive, define as 1 in this code
#AML negative, 0
#Given an unlabeled data point, you can make a prediction about what class 
#it belongs to based on the class of the data points in the training set that it is "closest" to. 

#The choice of distance metric depends on the data set you are working with. 
#We will be using Euclidean distance  as our distance metric.


#import necessary packages
import numpy
import sys
from random import shuffle


#function: check values of a list, if any of them is negative, replace it with zero
#input: two numpy matrix data points to be predicted and data points in the training set
#output: euclidean distance from each training data to 
####Euclidean###########
def cal_euclidean(pre_data,train_data):
    # calculate distance between each training data and all data to be predicted
    dst_out = []
    for train0 in train_data:
        #dst = distance
        dst = (pre_data - train0)**2
        dst = numpy.sum(dst,axis=1)
        dst = numpy.sqrt(dst)
        dst_out.append(dst)
    #transpose the output,so the output is in a form
    #n vectors (n = numbers of pre_data) of m long vectors (m = number of train_data)
    #each score is the distance between i(pre_data) and j (train_data)
    return numpy.array(dst_out).T

#function: sort list A, B based on A
def sort_based_2(list1,list2):
    list_two = zip(list1,list2)
    list1,list2 = zip(*sorted(list_two,key=lambda x: x[1]))
    return list1,list2

#function: check values of a list, if any of them is negative, replace it with zero
#input: two numpy matrix data points to be predicted and data points in the training set
#output: euclidean distance from each training data to 
def get_pred(dist0,train_label,k,p):
    index0 = range(0,len(dist0))
    [close_train,_] = sort_based_2(index0,dist0)
    #counting numbers of top postive near neighbors
    top_label = train_label[list(close_train[0:k])]
    pos_near_n = numpy.sum(top_label)
    #based on if more than p*k of near neighbors are positive
    #return 1 (postive) or 0
    if pos_near_n > p*k:
        return 1
    else:
        return 0

#fix
#function: check values of a list, if any of them is negative, replace it with zero
#input: a list
#output: a list of 
####KNN###########
def knn_predict(pre_data,train_data,train_label,k,p):
    #the output prediction
    pred_out = []
    #calculate Euclidean distance 
    dist_pre_to_train = cal_euclidean(pre_data, train_data)
    for dist0 in dist_pre_to_train:
        #get the top 
        result0 = get_pred(dist0,train_label,k,p)
        pred_out.append(result0)
    
    return pred_out



#function: calculate classifier performance based on true and predicted values
#input: list of predicted values (0,1) list of true values
#output: the following three parameters
#Sensitivity = TP/(TP+FN)
#Specificity = TN/(TN+FP)
#Accuracy = (TP+TN)/total
####Calculate performance###############
def cal_performance(list_pred, list_true):
    #initialize paramters 
    tp = 0
    tn = 0
    fp = 0
    fn = 0
    #verify the lengths of two lists
    if not len(list_pred) == len(list_true):
        print('Different length error')
    else:
        for i in range(0, len(list_pred)):
            result0 = (list_true[i],list_pred[i])
            if result0 == (0,0):
                tn +=1
            if result0 == (0,1):
                fp +=1
            if result0 == (1,0):
                fn +=1
            if result0 == (1,1):
                tp +=1
    #calculate the performance scores based on tp,tn,fp,fn
    sensi_i = float(tp)/(tp+fn)
    speci_i = float(tn)/(tn+fp)
    accu_i = float(tp+tn)/len(list_pred)
    #return performance
    return accu_i,sensi_i,speci_i


#fix
#function: check values of a list, if any of them is negative, replace it with zero
#input: a list
#output: a list with no negative scores
def get_data(file_name0,label0):
    data0 = []
    label_list = []
    for line0 in open(file_name0):
        line0 = line0.rstrip().split('\t')
        line0 = [float(w) for w in line0]
        data0.append(line0)
        label_list.append(label0)
    return numpy.array(data0).T, numpy.array(label_list)


#fix
#function: check values of a list, if any of them is negative, replace it with zero
#input: a list
#output: a list with no negative scores

def magic_round(x,n):
    return format(round(x,n),'.2f')

def write_output(file_name0,k,p,n,accu0,sensi0,speci0):
    fileout = open(file_name0,'w+')
    fileout.write('k: '+str(k)+'\n')
    fileout.write('p: '+magic_round(p,2)+'\n')
    fileout.write('n: '+str(n)+'\n')
    fileout.write('accuracy: '+magic_round(accu0,2)+'\n')
    fileout.write('sensitivity: '+magic_round(sensi0,2)+'\n')
    fileout.write('specificity: '+magic_round(speci0,2))
    fileout.close()

#fix
#function: check values of a list, if any of them is negative, replace it with zero
#input: a list
#output: a list with no negative scores
def knn_n_fold(pos_data,neg_data,pos_n_fold_set,neg_n_fold_set,k,p):
    accu0 = 0
    sensi0 = 0
    speci0 = 0
    #get the total data point number
    total_n = len(numpy.concatenate((pos_data,neg_data),axis=0))
    #get sets for pos and neg
    set_total_pos = set(range(0,len(pos_data)))
    set_total_neg = set(range(0,len(neg_data)))
    #n fold loop
    for i in range(0,len(pos_n_fold_set)):
        pos_train_index = list(set_total_pos - pos_n_fold_set[i])
        neg_train_index = list(set_total_neg - neg_n_fold_set[i])
        #generating training data
        train_data = numpy.concatenate([pos_data[pos_train_index],neg_data [neg_train_index]])
        #generate the labels (0,1) for the training 
        train_label = numpy.concatenate([numpy.ones(len(pos_train_index)),numpy.zeros(len(neg_train_index))])
        #generating validation data
        pred_data = numpy.concatenate([pos_data[list(pos_n_fold_set[i])],neg_data[list(neg_n_fold_set[i])]])
        #generate the true labels to calculate perforamnce
        list_true = numpy.concatenate([numpy.ones(len(pos_n_fold_set[i])),numpy.zeros(len(neg_n_fold_set[i]))])
        list_pred = knn_predict(pred_data,train_data,train_label,k,p)
        [accu_i,sensi_i,speci_i] = cal_performance(list_pred,list_true)
        #get validation number
        pred_n = len(pos_n_fold_set[i]) + len(neg_n_fold_set[i])
        #update performance scores
        accu0 += accu_i*pred_n
        sensi0 += sensi_i*pred_n
        speci0 += speci_i*pred_n
         
    #return weighted average values for three performance scores
    return accu0/total_n,sensi0/total_n,speci0/total_n

#fix
#function: check values of a list, if any of them is negative, replace it with zero
#input: a list
#output: a list with no negative scores
def get_n_fold(len0,n):
    list_out = []
    #calcualte div and mod in len0/n
    #to even out each group
    [div0,mod0] = divmod(len0,n)
    #starting index
    index0 = 0
    #n fold
    for _ in range(0,n):
        if mod0 > 0:
            #give the number one extra sample
            list_out.append(set(range(index0,index0+1+div0)))
            mod0 -=1 
            index0 += 1+div0
        else:
            list_out.append(set(range(index0,index0+div0)))
            index0 += div0
    #return a list of sets separting 0 to len0 to n parts with relatively even
    return list_out
            
#fix
#function: the main function of this project
#input: read in the input file with specification and sequence information
#output: perform Smith-Waterman alignment (global or local), output best matching scores
#        and all possible best alignments
def main():
    #initinize files and get I/O file names
    pos_file_name = sys.argv[1]
    neg_file_name = sys.argv[2]
    out_file_name = 'knn.out'
    k = int(sys.argv[3])
    p = float(sys.argv[4])
    #n fold validation
    n = int(sys.argv[5])
    #get all data 
    #all data and shuffle the data
    [pos_data,_] = get_data(pos_file_name, 1)
    shuffle(pos_data)
    #aml data and shuffle the data
    [neg_data,_] = get_data(neg_file_name, 0)
    shuffle(neg_data)
    #get n-fold validation
    pos_n_fold_set = get_n_fold(len(pos_data),n)
    neg_n_fold_set = get_n_fold(len(neg_data),n)
    #main prediction part
    [accu0,sensi0,speci0] = knn_n_fold(pos_data,neg_data,pos_n_fold_set,neg_n_fold_set,k,p)
    #output_performance into output file
    write_output(out_file_name,k,p,n,accu0,sensi0,speci0)
    

if __name__ == '__main__':
    main()        
        
        
        