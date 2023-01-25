from sympy import *
from sympy.stats import *
import numpy as np
import math
import time
from sympy.integrals.quadrature import gauss_laguerre


def lag_weights_roots(n):
    x = Symbol('x')
    roots = Poly(laguerre(n, x)).all_roots()
    x_i = [rt.evalf(20) for rt in roots]
    w_i = [(rt / ((n + 1) * laguerre(n + 1, rt)) ** 2).evalf(20) for rt in roots]
    return x_i, w_i


lagwr40 = lag_weights_roots(10)
print(lagwr40)
print(40)

beta = Symbol("beta")
theta = Symbol("theta")
x = Symbol("x")
u = Symbol("u")
lamb = Symbol("lamb")
alpha = 10 / beta
pxgut = density(Poisson("pxgut", theta * u))(x)
print(pxgut)
pt = density(Gamma("pt", alpha, beta))(theta)
print(pt)
print(pxgut * pt)
pxgu = integrate(pxgut * pt, (theta, 0, oo), conds='none')
qu = exp(-lamb * u) * lamb
print(pxgu * qu)
denomlog = integrate(qu * pxgu, (u, 0, oo), conds='none')
print(denomlog)
fnctn = (pxgu * qu * log(pxgu / denomlog, 2))
totalsollist = []
lenlist = len(lagwr40[0])
betalist = [1.5625]
for b in betalist:
    sollist = []
    for lambsub in range(41):
        L = 10 ** (lambsub / 40) / 5
        uintapprox40 = 0
        timeis = time.time()
        xsub = 0
        runpass = True
        xsave = 2000
        while xsub <= 2000 and runpass:
            funceva3 = (fnctn.subs(lamb, L).subs(x, xsub).subs(beta, b)) / exp(-u)
            valtoadd = 0
            for i in range(lenlist):
                valtoadd += N(lagwr40[1][i] * funceva3.subs(u, lagwr40[0][i]))
            change = valtoadd/(uintapprox40+valtoadd)
            usave = uintapprox40
            uintapprox40 += valtoadd
            if change < 1e-10 and change >0 and xsub >=50:
                xsave = xsub
                runpass = False
            xsub+=1
        print([uintapprox40,change,xsave, b, L])
        sollist.append([L, uintapprox40])
    totalsollist.append([b, sollist])
print(totalsollist)
