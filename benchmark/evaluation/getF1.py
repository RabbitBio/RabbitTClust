#!/usr/bin/env python3
from sklearn import metrics
import numpy as np
import sys
import pandas as pd

def classification_report_cvs(report, fileName):
    report_data = []
    lines = report.split('\n')
    for line in lines:
        report_data.append(line)
    #dataframe = pd.DataFrame.from_dict(report_data)
    #dataframe.to_csv(fileName+'.csv', index=False)

    line = lines[len(lines)-2] #for weighted precision, recall, and F1-score
    return line

def getF1(arr, fileName):
    prediction = arr[0,:]
    groundTruth = arr[1,:]
    report = metrics.classification_report(groundTruth, prediction, digits=3, zero_division=0)
    line = classification_report_cvs(report, fileName+'.F1')
    return line

if __name__ == "__main__":
    #fileList = ["rabbit_refseq", "rabbit_bacteria", "rabbit_half", "rabbit_sub", "mothur_refseq", "mothur_bacteria", "mothur_half", "mothur_sub", "mesh_half", "mesh_sub", "gclust_sub"];
    fileList = [sys.argv[1]]
    #fileList = ["rabbit_bacteria"]
    for file in fileList:
        originFile = np.loadtxt(file)
        result = getF1(originFile, file)
        print("result F1 of {} is: {} \n".format(file, result))
