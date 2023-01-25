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


def trial_receptor(x):
    Infoperk = {}
    info = {}
    #prove or discuss the ratios is all we need
    kmin = 0.01
    kstep = 0.5
    Kalist = [kmin*10**(i*kstep) for i in range(21)]
    Kblist = [kmin*10**(i*kstep) for i in range(21)]
    totalrunsneeded= len(Kblist)**2
    runcount = 0
    timestart = time.time()
    for kap in Kalist:
        for kbp in Kblist:


            kmin = min([kap,kbp,1])
            if kmin==0:
                kmin+=1
            Lmin = 100
            Lstepsize = 50
            Lstepnumber = 2
            Llist = [Lmin+Lstepsize*i for i in range(Lstepnumber)]
            #Use this to show rules still apply for weird distributions
            #Llist = [Lmin+Lstepsize*i for i in range(Lstepnumber) for j in range(i+1)]
            arrayofallpLRgL = []
            possibleLRlist = []
            pLlist = []
            acceptiblekeysforygx = []
            for Lv in Llist:
                Lvcount = 0
                for i in Llist:
                    if Lv == i:
                        Lvcount+=1
                pL = Lvcount/len(Llist)
                pLlist.append(pL)
            pLdict = dict(zip(Llist,pLlist))
            pLRgL = {}
            for Lv in Llist:
                class mRNA_Synth(Model):
                  def __init__(self, parameter_values = None):
                    # Initialize the model.
                    Model.__init__(self, name = "mass_action")

                    # Define parameters.
                    ka = Parameter(name = 'ka', expression = kap)
                    kb = Parameter(name='kb', expression=kbp)
                    ka2 = Parameter(name='ka2', expression=1)
                    kb2 = Parameter(name='kb2', expression=1)
                    self.add_parameter([ka,kb,ka2,kb2])

                    # Define molecular species.
                    L = Species(name = 'L', initial_value = Lv)
                    R = Species(name = 'R', initial_value = 100)
                    LR = Species(name='LR', initial_value=0)
                    LRs = Species(name='LRs', initial_value=0)
                    self.add_species([L,R,LR,LRs])

                    # Define reactions.
                    LandR_LR = Reaction(name = "LandR_LR", reactants = {L:1,R:1}, products = {LR:1},
                                  propensity_function = "L*R*ka")
                    LR_LandR = Reaction(name="LR_LandR", reactants={LR: 1}, products={L: 1, R: 1},
                                        propensity_function="LR*kb")
                    LR_LRs = Reaction(name="LR_LRs", reactants={LR: 1}, products={LRs: 1},
                                        propensity_function="LR*ka2")
                    LRs_LR = Reaction(name="LRs_LR", reactants={LRs: 1}, products={LRs: 1},
                                        propensity_function="LRs*kb2")
                    self.add_reaction([LandR_LR,LR_LandR,LR_LRs,LRs_LR])
                    self.timespan(numpy.linspace(0, 5*1/kmin,8))
                mRNA_Synth_model = mRNA_Synth()
                results = []
                n_results = mRNA_Synth_model.run(solver = gillespy2.solvers.NumPySSASolver, number_of_trajectories=10,debug = False)
                LRdist = []
                for i in n_results:
                    LRdist.append(i['LRs'][-1])
                    possibleLRlist.append(i['LRs'][-1])
                LRentries = list(set(LRdist))

                for i in LRentries:
                    LRcount = 0
                    for j in LRdist:
                        if i==j:
                            LRcount+=1
                    pLRgL[str([i,Lv])] = LRcount/len(LRdist)
                    acceptiblekeysforygx.append(str([i,Lv]))
            possibleLRlistsmall = list(set(possibleLRlist))
            pLRlistsmall = []
            #
            #Hey future andrew, play around with making L large enough so that you can basically assume L+R is fast enough
            #Hey future andrew, also remember L+R matters for the overall rate, try to predict the max occupancy and have your data bunch around that
            #
            pLR = {}
            for y in possibleLRlistsmall:
                py = 0
                for x in Llist:
                    if str([y,x]) in acceptiblekeysforygx:
                        print(str([y,x]))
                        py+=pLdict[x]*pLRgL[str([y,x])]
                pLR[y] = py
            print(pLR)
            ptest = 0
            for keys in pLR.keys():
                ptest+= pLR[keys]
            print(ptest)
            pLRlist = []
            pLRLdict = {}
            pptest = 0
            for keys in pLRgL.keys():
                pLRLdict[keys] = pLRgL[keys]*pLdict[eval(keys)[1]]
            for keys in pLRLdict.keys():
                pptest += pLRLdict[keys]
            print(pptest)
            print(pLRLdict)
            information = 0
            for y in possibleLRlistsmall:
                for x in Llist:
                    if str([y,x]) in acceptiblekeysforygx:
                        information += pLRLdict[str([y,x])]*math.log(pLRLdict[str([y,x])]/(pLdict[x]*pLR[y]),2)

            runcount+=1
            print("percent done: " + str(runcount/totalrunsneeded*100)+"%, time since start: " + str(time.time()-timestart) )
            Infoperk[str([kap,kbp])] = information
            info[str([kap, kbp])] = [pLRgL, pLR, pLdict, pLRLdict]
    print(Infoperk)
    print(info)
if __name__ == '__main__':
    timestart = time.time()
    trial_receptor(1)
    for repeat in range(1):
        with Pool(multiprocessing.cpu_count()) as p:
            #p.map(Noise_Maker,range(multiprocessing.cpu_count()))
            print(1)
    print(time.time()-timestart)