#!/usr/bin/env python
#BMI214 Project 2
#Binbin Chen
#bchen45

#Function
# Run K-means on the yeast microarray data. 
# We will cluster the genes based on their expression levels.

# Of the genes represented on the 79 microarrays, 121 were previously 
# characterized as ribosomal genes.Therefore, we might expect many
# of the ribosomal genes to be coordinately regulated -- they should have
# similar mRNA expression levels.

#Scanning yeast_gene_names.txt you will see that the ribosomes are the 
#last 121 genes in the file.

#you must number the clusters based on the order they are entered in 
#the centroid file, starting at 1. The first centroid in the file will
# be centroid 1, the 2nd will be centroid 2, etc.
# When centroids are generated randomly, they should be numbered 1 to k.

# the number of iterations completed should be
# written to <stdout> as shown below:
#iterations: 45

#import necessary packages
import numpy
import sys
import random
from collections import defaultdict
import matplotlib.pyplot as plt


#function: calculate euclidean distance between two sets of points (means/centroids and data)
#adopted from my KNN code
#input: two numpy matrix data points to be predicted and data points in the training set
#output: return the most nearby means from a given data point the same length of data point
####Euclidean###########
def compute_closest_means(means,data):
    out = []
    for m in means:
        dst = (data - m)**2
        dst = numpy.sum(dst,axis=1)
        dst = numpy.sqrt(dst)
        out.append(dst)
    #return the most nearby means
    return numpy.argmin(numpy.array(out).T, axis=1)

#fix
#function: check values of a list, if any of them is negative, replace it with zero
#input: a list
#output: a list with no negative scores
def compute_new_means(closest_means, data):
    groups = defaultdict(list)
    for i in range(len(closest_means)):
        #fix mer
        groups[closest_means[i]].append(data[i])
    means = []
    for k,v in sorted(groups.iteritems(),key=lambda x: x[0]):
        means.append(numpy.average(v,axis=0))
    return numpy.array(means)

#fix
#function: check values of a list, if any of them is negative, replace it with zero
#input: a list
#output: a list with no negative scores
def same_means(samp_one, samp_two):
    match_dict = {}
    for x,y in zip(list(samp_one),list(samp_two)):
        if x in match_dict:
            if y != match_dict[x]:
                return False
        else:
            match_dict[x] = y
    return True

#fix
#function: check values of a list, if any of them is negative, replace it with zero
#input: a list
#output: a list with no negative scores
def kmeans(k,data,means,max_it):
    for i in range(1,max_it+1):
        #i is the iteration number 
        new_closest_means = compute_closest_means(means, data)
        if (i>1) and (same_means(new_closest_means,closest_means)):
            #this is virtually a break 
            return closest_means, i
        else:
            closest_means = new_closest_means
        means = compute_new_means(closest_means, data)
    return closest_means, i


#fix
#function: check values of a list, if any of them is negative, replace it with zero
#input: a list
#output: a list with no negative scores
def write_output(file_name0,k,p,n,accu0,sensi0,speci0):
    fileout = open(file_name0,'w+')
    fileout.write('k: '+str(k)+'\n')
    fileout.write('p: '+str(round(p,2))+'\n')
    fileout.write('n: '+str(n)+'\n')
    fileout.write('accuracy: '+str(round(accu0,2))+'\n')
    fileout.write('sensitivity: '+str(round(sensi0,2))+'\n')
    fileout.write('specificity: '+str(round(speci0,2)))
    fileout.close()


#fix
#function: the main function of this project
#input: read in the input file with specification and sequence information
#output: perform Smith-Waterman alignment (global or local), output best matching scores
def generate_centro(data0,k):
    data0 = data0.T
    data_centro = []
    for i in range(0,k):
        centro0 = []
        for j in range(0,len(data0)):
            min0 = numpy.min(data0[j])
            max0 = numpy.max(data0[j])
            diff0 = max0 - min0
            value0 = min0 + random.random()*diff0
            centro0.append(value0)
        data_centro.append(centro0)
    return numpy.array(data_centro)

#fix
#function: check values of a list, if any of them is negative, replace it with zero
#input: a list
#output: a list with no negative scores
def get_data(file_name0):
    data0 = []
    for line0 in open(file_name0):
        line0 = line0.rstrip().split('\t')
        line0 = [float(w) for w in line0]
        data0.append(line0)
    return numpy.array(data0)


#fix
#function: the main function of this project
#input: read in the input file with specification and sequence information
#output: perform Smith-Waterman alignment (global or local), output best matching scores
#        and all possible best alignments
def main():
    #initinize files and get I/O file names
    '''
    k = sys.argv[1]
    data_file_name = sys.argv[2]
    out_file_name = 'knn.out'
    #get data from data_file_name
    data0 = get_data(path0+data_file_name)
    # maximum iteration
    max_it = sys.argv[3]
    # if available get centroid files and load data
    # if not, randomly generating centroids
    if len(sys.argv) > 4:
        centro_file_name = sys.argv[4]
        data_centro = get_data(path0+centro_file_name)
        data_centro = data_centro[:k]
    else:
        data_centro = random.sample(list(data0),k)
    '''
        
    #for debugging on local only
    #data_file_name = '/testData/testdata.dat'
    #centro_file_name = '/testData/testdata_centroids.dat'
    data_file_name = 'yeastData/yeast.dat'
    centro_file_name = 'yeastData/yeast_centroids.txt'
    data0 = get_data(path0+data_file_name)
    data_centro = get_data(path0+centro_file_name)
    k = 2
    #data_centro = generate_centro(data0,k)
    data_centro = data_centro[:k]
    data_centro = []
    data_centro.append(data0[0])
    data_centro.append(data0[2466])
    data_centro = numpy.array(data_centro)
    max_it  = 100
    data_centro = random.sample(list(data0),k)
    #neg_file_name = 'AML.dat'
    #out_file_name = 'knn.out'
    #k = 10
    #p = 0.5
    #n = 4
    #get all data 
    #plotting
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #ax.plot(listx,listy,'ro-')
    #data_x = data0.T[0]
    #data_y = data0.T[1]
    #ax.plot(data_x,data_y,'*r')
    #fig.savefig(path0+'plot_test.png')
    
    
    [list_out, iter0] = kmeans(k,data0,data_centro,max_it)
    print(len(list_out))
    print(list_out[-121:],list_out[1510])
    print('iteration')
    print(iter0)
    

    
    


#the following line is for debugging in my local machine only
path0 = '/Users/binbineow2/Dropbox/Russ/python_code/Common_Work_Station/BMI214/kmeans/'
#input_name = path0+'alignment_example2.input'
#output_name = path0+'alignment_example2.output_user'
#path0 = ''

if __name__ == '__main__':
    main()        
        
        
        