import numpy as np
import math
import scipy.special
import time
import matplotlib.pyplot as plt

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
    return [val,(t1-t2)/t1]
totalstorelist = []
#betalist = [0.004,0.036,0.1,0.196,0.324,0.484,0.676,0.9,1.156,1.444,1.764]
betalist = [20]
tstart = time.time()
mcacc = 1000
mcrec = 1000
count = 0
for b in betalist:
    sollist = []
    for lambsub in range(1):
        timestart = time.time()
        L = 10 ** (lambsub / 40)/5
        scale = b  # mean=4, std=2*sqrt(2)
        shape = 10/scale
        sumnp = [0 for i in range(mcrec)]
        thetalist = np.random.gamma(shape, scale, size=mcacc)
        # for error analysis
        #plt.hist(thetalist, bins=100)
        #plt.show()
        #thetalist = np.ones(mcacc)
        # for error analysis
        for i in range(len(sumnp)):

            totalval = 0
            for x in range(2000):

                u1list = np.random.exponential(scale=1 / L, size=mcacc)
                pxgut1 = np.exp(x*np.log(np.multiply(thetalist,u1list))-np.multiply(thetalist,u1list)-scipy.special.loggamma(x+1))
                denomlog = np.average(pxgut1)
                u2list = np.random.exponential(scale=1 / L, size=mcacc)
                pxgut2 = np.exp(x*np.log(np.multiply(thetalist,u2list))-np.multiply(thetalist,u2list)-scipy.special.loggamma(x+1))
                value = pxgut2*np.log2(pxgut2/denomlog)
                value[np.isnan(value)] = 0
                value[np.isinf(value)] = 0
                #value = np.average(value)
                totalval += value
                #if x % 1000 == 0:
                    #print([b,lambsub,np.average(value),np.average(totalval)])
            sumnp[i] = np.average(totalval)
            print([b, lambsub,i,np.average(totalval)] )
        npsunp = np.array(sumnp)
        errortotalval = np.std(npsunp)
        averagetotal = np.average(npsunp)
        print([L,averagetotal,errortotalval,-1*(timestart - time.time())/60*(41*11-count)])
        count+=1
        sollist.append([L,averagetotal,errortotalval])
    totalstorelist.append([b,sollist])
print(totalstorelist)
print(time.time()-tstart)
input()
for i in totalstorelist:
    print(i)
    print("i")
    for j in i[1]:
        print([i[0],j[0],j[1]])
        print("j")
    CV = (i[0]/10)**0.5
    name = "MI" + str(CV) + ".txt"
    with open(name, "a+") as file:
        for j in i[1]:
            file.write(str(j[0]) + " " + str(j[1]) + " " + str(j[2]) + "\n")
            file.flush()




