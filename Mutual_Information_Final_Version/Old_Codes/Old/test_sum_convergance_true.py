import numpy as np
import math
import scipy.special
import time

def logfac(x):
    lf = 0
    for i in range(1,x+1):
        lf+= np.log(i)
    return lf
def singleseqval(theta,lamb,x):
    return theta ** x * lamb * (theta + lamb) ** (-x - 2) * (-theta - x * theta + (theta + lamb) * (
                -np.log(lamb) + np.log(theta + lamb) - np.log(scipy.special.gamma(1 + x)) + x * scipy.special.digamma(
            1 + x))) / np.log(2)
def singleseqval2(theta,lamb,x):
    return (theta/(theta+lamb)) ** x * lamb * (theta + lamb) ** ( - 2) * (-theta - x * theta + (theta + lamb) * (
                -np.log(lamb) + np.log(theta + lamb) - logfac(x) + x * scipy.special.digamma(
            1 + x))) / np.log(2)
def getvalue(theta,lamb,maxX):
    running = True
    val = 0
    for x in range(maxX):
        val += singleseqval2(theta,lamb,x)
        #test = val+singleseqval(theta,lamb,x+1)/(1-singleseqval(theta,lamb,x+1)/singleseqval(theta,lamb,x))
    return val
def listgetval(thetalist,lamb,maxX):
    val = thetalist-thetalist
    for x in range(maxX):
        val += (thetalist/(thetalist+lamb)) ** x * lamb * (thetalist + lamb) ** ( - 2) * (-thetalist - x * thetalist + (thetalist + lamb) * (
                -np.log(lamb) + np.log(thetalist + lamb) - logfac(x) + x * scipy.special.digamma(
            1 + x))) / np.log(2)
        if x % (maxX-25) ==0:
            valsave = val.copy()
    t1 = 0
    for i in val:
        t1 += i
    t1 /= len(val)
    t2 = 0
    for i in valsave:
        t2 += i
    t2 /= len(valsave)
    print("error indicator: " + str((t1-t2)/t1))
    return val
totalstorelist = []
betalist = [0.004,0.036,0.1,0.196,0.324,0.484,0.676,0.9,1.156,1.444,1.764]
for b in betalist:
    sollist = []
    for lambsub in range(41):
        L = 10 ** (lambsub / 40)/5
        timestart = time.time()
        scale = b  # mean=4, std=2*sqrt(2)
        shape = 10/scale
        s = np.array([10.0])
        thetalist = s
        pthetalist = s
        test = 0
        value = listgetval(thetalist,L,2000)
        for i in thetalist:
            test+= i
        t= 0
        for i in value:
            t+= i
        t/=len(value)
        print([L,t])
        sollist.append([L,t])
    totalstorelist.append([b,sollist])
print(totalstorelist)
for i in totalstorelist:
    print(i)
    print("i")
    for j in i[1]:
        print([i[0],j[0],j[1]])
        print("j")
    CV = (i[0]/10)**0.5
    name = "TrueMI" + str(CV) + ".txt"
    with open(name, "a+") as file:
        for j in i[1]:
            file.write(str(j[0]) + " " + str(j[1]) + "\n")
            file.flush()

