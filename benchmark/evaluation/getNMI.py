#!/usr/bin/env python3
from sklearn import metrics
import numpy as np
import sys

def getNMI(arr):
    a = arr[0,:]
    b = arr[1,:]
    result = metrics.normalized_mutual_info_score(a, b)
    return result

if __name__ == "__main__":
    #fileList = ["rabbit_refseq", "rabbit_bacteria", "rabbit_half", "rabbit_sub", "mothur_refseq", "mothur_bacteria", "mothur_half", "mothur_sub", "mesh_half", "mesh_sub", "gclust_sub"];
    fileList = [sys.argv[1]]
    #fileList = ["rabbit_bacteria"]
    for file in fileList:
        originFile = np.loadtxt(file)
        result = getNMI(originFile)
        print("result NMI of {} is: {} \n".format(file, result))
