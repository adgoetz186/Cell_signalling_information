import numpy
import testwhat
import gillespy2
from gillespy2 import Model, Species, Parameter, Reaction
from gillespy2.solvers.numpy.basic_ode_solver import BasicODESolver
import matplotlib.pyplot as plt
import time
import math
import os
from multiprocessing import Pool
import multiprocessing
global mRNAProtein

with open("nmkr", "r") as runDataFile:
    mRNAProtein = eval(runDataFile.read())
    print(testwhat.choice(mRNAProtein))
def initalize(x):
    mRNAandProteinCombo = testwhat.choice(mRNAProtein)
    mRNA = mRNAandProteinCombo[0]
    protein =

if __name__ == '__main__':
    #os.chdir('/blue/pdixit/agoetz/Raw_Output/')
    print(mRNAProtein[0])
    mRNA = [sublist[0] for sublist in mRNAProtein]
    protein = [sublist[1] for sublist in mRNAProtein]
    kmprod = math.log(2) / ((8.25) * (60 ** 2)) * 10.44
    kmdeg = math.log(2) / ((8.25) * (60 ** 2))
    knprod = 21.58 / (60 ** 2)
    kndeg = math.log(2) / ((22.23) * (60 ** 2))
    mRNAmean = 0
    for i in mRNA:
        mRNAmean += i
    mRNAmean /= len(mRNA)
    proteinmean = 0
    for i in protein:
        proteinmean += i
    proteinmean /= len(protein)
    print(len(mRNAProtein))
    print("predicted mRNA mean is: " + str(kmprod / kmdeg) + ", actual mRNA mean is: " + str(mRNAmean))
    print("predicted protein mean is: " + str(kmprod / kmdeg * knprod / kndeg) + ", actual protein mean is: " + str(
        proteinmean))
    timestart = time.time()
    with Pool(multiprocessing.cpu_count()) as p:
        output = p.map(initalize,range(multiprocessing.cpu_count()))
