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


# the number of iterations completed should be
# written to <stdout> as shown below:
#iterations: 45

#import necessary packages
import numpy
import sys
import random
from collections import defaultdict


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


#function: update the new centroids
#input: a list of distance and data points
#output:reutnr the new centroids/means matrix
def compute_new_means(closest_means, data):
    groups = defaultdict(list)
    for i in range(len(closest_means)):
        groups[closest_means[i]].append(data[i])
    means = []
    for k,v in sorted(groups.iteritems(),key=lambda x: x[0]):
        means.append(numpy.average(v,axis=0))
    return numpy.array(means)


#function: check two list of clustering to see if they are the same
#input: two lists of cluster assigment
#output: whether if any assignment is updates
def same_means(samp_one, samp_two):
    match_dict = {}
    for x,y in zip(list(samp_one),list(samp_two)):
        if x in match_dict:
            if y != match_dict[x]:
                return False
        else:
            match_dict[x] = y
    return True


#function: main kmeans function to assign cluster until it convers or reach the maximum iteration
#input: data, centroids, and paramters k and max_int
#output: assignments and iteration when the kmeans ended
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



#function: ouptut performance to knn.out
#input: assignment of each points to clusters
#output: write values into the kmeans.out
def write_output(file_name0,list0):
    fileout = open(file_name0,'w+')
    index0 = 1
    for cluster0 in list0:
        fileout.write(str(index0)+'\t'+str(cluster0+1)+'\n')
        index0 +=1

#function: get data from file_name0 based on label0 to dtermine it's ALL or AML
#input: file_name string 
#output: numpy array of all data
def get_data(file_name0):
    data0 = []
    for line0 in open(file_name0):
        line0 = line0.rstrip().split('\t')
        line0 = [float(w) for w in line0]
        data0.append(line0)
    return numpy.array(data0)



#function: the main function of this project
#input: read in the input file with k value
#output: clusters of points based on centroid provided or generated
def main():
    #initinize files and get I/O file names
    k = int(sys.argv[1])
    data_file_name = sys.argv[2]
    out_file_name = 'kmeans.out'
    #get data from data_file_name
    data0 = get_data(data_file_name)
    # maximum iteration
    max_it = int(sys.argv[3])
    # if available get centroid files and load data
    # if not, randomly generating centroids
    if len(sys.argv) > 4:
        centro_file_name = sys.argv[4]
        data_centro = get_data(centro_file_name)
        data_centro = data_centro[:k]
    else:
        data_centro = random.sample(list(data0),k)

    [list_out, iter0] = kmeans(k,data0,data_centro,max_it)
    #print iteration when it converges or stops
    print('iterations: '+str(iter0))
    #write output file knn.out
    write_output(out_file_name,list_out)



if __name__ == '__main__':
    main()        
        
        
        