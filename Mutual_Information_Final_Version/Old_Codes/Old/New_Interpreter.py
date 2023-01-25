import numpy
numpy.max_line_width=numpy.inf
import gillespy2
from gillespy2 import Model, Species, Parameter, Reaction
from gillespy2.solvers.numpy.basic_ode_solver import BasicODESolver
import matplotlib.pyplot as plt
import time
import math
import os
from gekko import GEKKO
from multiprocessing import Pool
import multiprocessing
import json

def stringDictToDict(stringdict):
    outputarraykeys = stringdict.replace("'", '"').split('"')
    outputarrayvalues = stringdict.replace(".", "").replace("(", ')').split(')')
    keylist = []
    valuelist = []
    for i in range(len(outputarraykeys)):
        if i % 2 == 1:
            keylist.append(outputarraykeys[i])
    for i in range(len(outputarrayvalues)):
        if i % 2 == 1:
            valuelist.append(eval(outputarrayvalues[i]))
    rundictionary = dict(zip(keylist, valuelist))
    return rundictionary

def runliststoaveragelist(rlists):
    averagelist = [0 for i in range(len(rlists[0]))]
    for i in rlists:
        for j in range(len(i)):
            averagelist[j] += i[j]
    for i in range(len(averagelist)):
        averagelist[i] /= len(rlists)
    print(averagelist)
    return averagelist




filename = "Raw_Output2/Fold_is_2____a_is_-2.31____b_is_-2.5.txt"
if os.path.exists(filename):
    with open(filename, "r") as runDataFile:
        alllines = runDataFile.readlines()
    runlist = []
    for i in alllines:
        runlist.append(stringDictToDict(i))
    timelist = runlist[0]['time']
    phosOut = runlist[0]['sLR2']
    phosOutfirstmoment = phosOut.copy()
    phosOutsecondmoment = phosOut.copy()
    for i in range(len(phosOutsecondmoment)):
        phosOutsecondmoment[i] *= phosOutsecondmoment[i]
    for i in range(1,len(runlist)):
        for j in range(len(phosOutfirstmoment)):
            phosOutfirstmoment[j] += runlist[i]['sLR2'][j]
            phosOutsecondmoment[j] += runlist[i]['sLR2'][j]**2
    for i in range(len(phosOutfirstmoment)):
        phosOutfirstmoment[i]/=len(runlist)
        phosOutsecondmoment[i] /= len(runlist)
    std = [0 for i in range(len(phosOutfirstmoment))]
    for i in range(len(phosOutfirstmoment)):
        std[i] = phosOutsecondmoment[i] - phosOutfirstmoment[i]**2
    plt.plot(timelist,phosOutfirstmoment,'g')
    plt.show()
    plt.plot(timelist,phosOutsecondmoment,'r')
    plt.show()
    plt.plot(timelist, std, 'b')
    plt.show()







else:
    print("File not found")