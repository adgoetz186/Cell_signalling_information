import numpy as np
import math
import scipy.special
import time
#betalist = [0.004,0.036,0.1,0.196,0.324,0.484,0.676,0.9,1.156,1.444,1.764]
betalist = [0.004,0.036,0.1,0.196,0.324,0.484,0.676,0.9,1.156,1.444,1.764]
betalist = [20]
klist = [0.5,0,-0.5]
bc = 0
tstart = time.time()
mcacc = 10000
mcrec = 10
count = 0
fulldatalist = []
for mean in range(3):
    meanval = 10 ** klist[mean]
    datalist = []
    for b in betalist:
        sollist = []
        for inputscale in range(9):
            timestart = time.time()
            inputscalevalue = 10 ** (inputscale / 8)/5
            scale = b  # mean=4, std=2*sqrt(2)
            shape = 10/scale
            sumnp = [0 for i in range(mcrec)]
            shapeinput = meanval ** 2 * inputscalevalue ** 2
            scaleinput = meanval / shapeinput

            thetalist = np.random.gamma(shape, scale, size=mcacc)
            # for error analysis
            thetalistval = np.random.gamma(shape, scale, size=mcacc)
            thetalist = np.ones(mcacc) * thetalistval[0]
            for j in range(1,len(thetalistval)):
                addlist =  np.ones(mcacc) * thetalistval[j]
                thetalist = np.concatenate((thetalist,addlist),axis = None)
            # for error analysis
            totalval = 0
            for x in range(2000):
                u1list = np.random.gamma(shape=shapeinput,scale = scaleinput, size=mcacc*mcacc)
                pxgut1 = np.exp(x*np.log(np.multiply(thetalist,u1list))-np.multiply(thetalist,u1list)-scipy.special.loggamma(x+1))
                denomlog = np.average(pxgut1)
                u2list = np.random.gamma(shape=shapeinput,scale = scaleinput, size=mcacc*mcacc)
                pxgut2 = np.exp(x*np.log(np.multiply(thetalist,u2list))-np.multiply(thetalist,u2list)-scipy.special.loggamma(x+1))
                value = pxgut2*np.log2(pxgut2/denomlog)
                value[np.isnan(value)] = 0
                value[np.isinf(value)] = 0
                #value = np.average(value)
                totalval += value
                if x % 100 == 0:
                    print([b,np.average(value),np.average(totalval)])
            npsunp = np.array(totalval)
            errortotalval = np.std(npsunp)
            averagetotal = np.average(npsunp)
            print([mean,b,inputscale,averagetotal,errortotalval,-1*(timestart - time.time())/60*(9*3*11-count)])
            count+=1
            sollist.append([inputscale,averagetotal,errortotalval,shapeinput,scaleinput])
        datalist.append([b,sollist])
    fulldatalist.append([meanval,datalist])

print(fulldatalist)
print(time.time()-tstart)





