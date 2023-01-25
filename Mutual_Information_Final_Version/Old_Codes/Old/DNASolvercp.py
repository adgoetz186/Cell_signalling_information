from sympy import *
from sympy.stats import *
import numpy as np
import math
import time
from sympy.integrals.quadrature import gauss_laguerre

def lag_weights_roots(n):
    x = Symbol('x')
    roots = Poly(laguerre(n, x)).all_roots()
    x_i = [rt.evalf(50) for rt in roots]
    w_i = [(rt/((n+1)*laguerre(n+1, rt))**2).evalf(50) for rt in roots]
    return x_i, w_i








lagwr40 = lag_weights_roots(4)
print(lagwr40)
print(40)
lagwr40 = ([0.035700394308888385122, 0.18816228315869851600, 0.46269428131457645357, 0.85977296397293492226, 1.3800108205273371865, 2.0242091359228267334, 2.7933693535068164577, 3.6887026779082702096, 4.7116411465549726936, 5.8638508783437181143, 7.1472479081022882507, 8.5640170175861637627, 10.116634048451939407, 11.807892294004584843, 13.640933712537087228, 15.619285893339073837, 17.746905950095663043, 20.028232834574890530, 22.468249983498418351, 25.072560772426203794, 27.847480009168862721, 30.800145739445462701, 33.938657084913719609, 37.272245880476004328, 40.811492823886920466, 44.568603175334462707, 48.557763533059992281, 52.795611187216932969, 57.301863323393627495, 62.100179072775111612, 67.219370927126998799, 72.695158847612462118, 78.572802911571309281, 84.911231135704984543, 91.789874671236376992, 99.320808717446808250, 107.67244063938827252, 117.12230951269068881, 128.20184198825565119, 142.28004446915999789], [0.088412106190342440940, 0.17681473909572229560, 0.21136311701596243103, 0.19408119531860179966, 0.14643428242412511441, 0.093326798435770880507, 0.050932204361044237026, 0.023976193015684841840, 0.0097746252467144596189, 0.0034579399930184868613, 0.0010622468938968719350, 0.00028327168532432471583, 0.000065509405003246292798, 0.000013116069073267784125, 2.2684528787793650545e-6, 3.3796264822006792108e-7, 4.3228213222820885689e-8, 4.7284937709907793279e-9, 4.4031741042328488129e-10, 3.4724414848038224856e-11, 2.3053815449168221616e-12, 1.2797725976766356072e-13, 5.8941771723511529447e-15, 2.2322175799045774184e-16, 6.8803364842843023409e-18, 1.7056037368180867485e-19, 3.3537119406661829355e-21, 5.1461995601366791408e-23, 6.0447625115876632890e-25, 5.3105847773213307528e-27, 3.3925280532805218961e-29, 1.5217354931814569975e-31, 4.5852916145026869176e-34, 8.7621586574862485610e-37, 9.8274157251479333061e-40, 5.8011520191697791085e-43, 1.5309086846066868536e-46, 1.3819863056493280997e-50, 2.5666336050123721838e-55, 2.7003609402170336406e-61])




totalsollist = []
for betasub in range(3):
    b = 10**(betasub-1)
    sollist = []
    for lambd in range(200):
        L = (lambd+1)/100
        theta = Symbol("theta")
        x = Symbol("x")
        u = Symbol("u")
        beta = b
        lamb = L
        alpha = 1 / beta
        pxgut = density(Poisson("pxgut", theta * u))(x)
        pt = density(Gamma("pt", alpha, beta))(theta)
        pxgu = integrate(pxgut * pt, (theta, 0, oo), conds='none')
        qu = exp(-lamb * u) * lamb
        denomlog = integrate(qu * pxgu, (u, 0, oo), conds='none')
        fnctn = (pxgu * log(pxgu / denomlog, 2))
        #Note the above is the function without the qu (lamb*exp(-lamb * u))
        uintapprox40 = 0
        for xsub in range(10000):
            funceva3 = fnctn.subs(lamb, L).subs(x, xsub).subs(beta, b)
            valtoaddmont = 0
            monttime = time.time()
            explist = np.random.exponential(1/L, 50)
            for i in explist:
                valtoaddmont += N(funceva3.subs(u,np.log(i*1/L ))*-1/L)
            valtoaddmont /= len(explist)
            monttime = time.time() - monttime
            valtoadd = 0
            for i in range(len(lagwr40[0])):
                valtoadd += N(lagwr40[1][i] * funceva3.subs(u, lagwr40[0][i]))
                if funceva3.subs(u, lagwr40[0][i]) < 0:
                    print([funceva3.subs(u, lagwr40[0][i]),N(funceva3.subs(u, lagwr40[0][i])),lagwr40[0][i]])
            input([valtoaddmont,monttime,valtoadd])
            usave = uintapprox40
            uintapprox40 += valtoaddmont
            if usave != 0:
                print((uintapprox40-usave)/usave)
                print(uintapprox40)
                if (uintapprox40-usave)/usave <1e-7 and (uintapprox40-usave)/usave>0and xsub >20:
                    break
        print([L,b])
        sollist.append([L,uintapprox40])
    totalsollist.append([b,sollist])
print(totalsollist)
