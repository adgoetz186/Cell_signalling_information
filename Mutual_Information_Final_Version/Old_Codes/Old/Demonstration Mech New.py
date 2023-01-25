import numpy as np
import math
import scipy.special as sc
import time
import matplotlib.pyplot as plt


totalstorelist = []
betalist = [0.004,0.036,0.1,0.196,0.324,0.484,0.676,0.9,1.156,1.444,1.764]
#betalist = [20]
tstart = time.time()
mcacc = 100000
mcrec = 100
count = 0



meanval = 10 ** (1/2)
inputscalevalue = 10 ** (0) / 5
shapeinput = meanval ** 2 * inputscalevalue ** 2
scaleinput = meanval / shapeinput


tlist = []




thetalist = np.ones(mcacc) * 10
totalval1 = 0
for x in range(2000):
    u1list = np.ones(mcacc)*10**(1/2)
    u2list = np.ones(mcacc)*10**(1/2)
    pxguttest = np.exp(x * np.log(np.multiply(thetalist, u1list)) - np.multiply(thetalist, u1list) - sc.loggamma(x + 1))
    logpxgut1 = x*np.log(np.multiply(thetalist,u1list))-np.multiply(thetalist,u1list)-sc.loggamma(x+1)
    logdenomlog = sc.logsumexp(logpxgut1) - np.log(mcacc)

    logpxgut2 = x*np.log(np.multiply(thetalist,u2list))-np.multiply(thetalist,u2list)-sc.loggamma(x+1)
    pxgut2 = np.exp(logpxgut2)

    pxgut2test = np.exp(x * np.log(np.multiply(thetalist, u2list)) - np.multiply(thetalist, u2list) - sc.loggamma(x + 1))



    value = pxgut2*logpxgut2 - pxgut2*logdenomlog

    totalval1 += value


print(totalval1)





