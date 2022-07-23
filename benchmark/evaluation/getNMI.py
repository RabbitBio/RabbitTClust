#!/usr/bin/env python3
from sklearn import metrics
import numpy as np
import sys
import pandas as pd

def classification_report_cvs(report):
    report_data = []
    lines = report.split('\n')
    #for line in lines:
    #    report_data.append(line)
    #dataframe = pd.DataFrame.from_dict(report_data)
    #dataframe.to_csv(fileName+'.csv', index=False)
    line = lines[len(lines)-2] #for weighted precision, recall, and F1-score
    return line

def getF1(arr):
    prediction = arr[0,:]
    groundTruth = arr[1,:]
    report = metrics.classification_report(groundTruth, prediction, digits=3, zero_division=0)
    line = classification_report_cvs(report)
    return line

def getNMI(arr):
    a = arr[0,:]
    b = arr[1,:]
    result = metrics.normalized_mutual_info_score(a, b)
    return result

if __name__ == "__main__":
    fileList = [sys.argv[1]]
    for file in fileList:
        originFile = np.loadtxt(file)
        F1 = getF1(originFile)
        NMI = getNMI(originFile)
        print("result F1 of {} is: {} \n".format(file, F1))
        print("result NMI of {} is: {} \n".format(file, NMI))
