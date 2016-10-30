import numpy
import random
from collections import defaultdict

test_data = numpy.random.randint(low=0,high=10,size=(100,50))

def compute_closest_means(means,data):
    out = []
    for m in means:
        dst = (data - m)**2
        dst = numpy.sum(dst,axis=1)
        dst = numpy.sqrt(dst)
        out.append(dst)
    return numpy.argmin(numpy.array(out).T, axis=1)

def compute_new_means(closest_means, data):
    groups = defaultdict(list)
    for i in range(len(closest_means)):
        groups[closest_means[i]].append(data[i])
    means = []
    for k,v in groups.items():
        means.append(numpy.average(v,axis=0))
    return numpy.array(means)

def same_means(samp_one, samp_two):
    match_dict = {}
    for x,y in zip(list(samp_one),list(samp_two)):
        if x in match_dict:
            if y != match_dict[x]:
                return False
        else:
            match_dict[x] = y
    return True

def kmeans(k,data,num_iter=1000):
    if k > data.shape[0]:
        print("Can't find {} means for {} data points".format(k,data.shape[0]))
        return None
    means = random.sample(list(data),k) #sample doesn't work with numpy arrays
    for i in range(num_iter):
        new_closest_means = compute_closest_means(means, data)
        if (i>0) and (same_means(new_closest_means,closest_means)):
            return closest_means
        else:
            closest_means = new_closest_means
        means = compute_new_means(closest_means, data)
    return closest_means

test_data = numpy.array([[-5, -3], [-3, -4], [-4,-5], [-6,3], [-5,5], [-3,6], [3,7], [4,5], [5,5]])

kmeans(3, test_data)
