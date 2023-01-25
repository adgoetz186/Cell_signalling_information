from sympy import *
from sympy.stats import *
import numpy as np
import math
import time
import scipy as sp
import scipy.special as sc
from decimal import *
def logbeta(a,b):
    return sc.loggamma(a)+sc.loggamma(b)-sc.loggamma(a+b)
lgvalues = np.polynomial.laguerre.laggauss(40)
datalist = []
betalist = [0.004,0.036,0.1,0.196,0.324,0.484,0.676,0.9,1.156,1.444,1.764]

bc = 0
count = 0
for b in betalist:
    k = 10/b
    betaentry = []
    for lambsub in range(41):
        timestart = time.time()
        L = 10 ** (lambsub / 40) / 5
        mcrecurs = 100
        mcacc = 1000
        sumnp = [0 for i in range(mcrecurs)]
        for i in range(len(sumnp)):
            for x in range(2000):
                u1list = np.random.exponential(scale=1 / L, size=mcacc)
                if x == 0:
                    pxgu1 = np.exp(x*np.log(u1list)-x*np.log(u1list + 1 / b)+k*np.log(1 / b * 1 / (u1list + 1 / b)))
                    denomlog = np.average(pxgu1)
                    u2list = np.random.exponential(scale=1 / L, size=mcacc)
                    pxgu2 = np.exp(x * np.log(u2list) - x * np.log(u2list + 1 / b) + k*np.log(1 / b * 1 / (u2list + 1 / b)))
                    integralval = pxgu2*np.log2(pxgu2/denomlog)
                    integralval[np.isnan(integralval)] = 0
                    integralval = np.average(integralval)
                else:
                    pxgu1 = np.exp(-1 * np.log(x) -  logbeta(k, x)+k*np.log(1 / b * 1 / (u1list + 1 / b))+ x*np.log(u1list / (u1list + 1 / b)))
                    denomlog = np.average(pxgu1)
                    u2list = np.random.exponential(scale=1 / L, size=mcacc)
                    pxgu2 = np.exp(-1 * np.log(x) -  logbeta(k, x)+k*np.log(1 / b * 1 / (u2list + 1 / b))+x*np.log(u2list/(u2list + 1 / b)))
                    integralval = pxgu2 * np.log2(pxgu2 / denomlog)
                    integralval[np.isnan(integralval)] = 0
                    integralval[np.isinf(integralval)] = 0
                    integralval = np.average(integralval)
                sumnp[i] += integralval
        npsunp = np.array(sumnp)
        betaentry.append([L,np.average(npsunp),np.std(npsunp)])
        count +=1
        print([b,L,lambsub+1,np.average(npsunp),np.std(npsunp),np.std(npsunp)/np.average(npsunp),(time.time()-timestart)/60,(time.time()-timestart)/60*((41)*(11)-count),41*11-count])
    datalist.append([b,betaentry])
    bc +=1
print(datalist)


