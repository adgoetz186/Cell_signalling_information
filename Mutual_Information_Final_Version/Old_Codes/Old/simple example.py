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

plt.hist([17.0],bins=50)
plt.xlabel("single run protein level at about 12 hours")
plt.ylabel("Number of times this value was recorded during 1 run")
plt.show()
def trial_receptor(x):
    Infoperk = {}
    info = {}
    #prove or discuss the ratios is all we need
    kstep = 0.5
    Kalist = [10**(3+i*kstep) for i in range(16)]
    Kblist = [10**(-5+i*kstep) for i in range(10)]
    Kalist = [10**7.5]
    Kblist = [10**(-3)]
    totalrunsneeded= len(Kblist)*len(Kalist)
    runcount = 0
    timestart = time.time()
    for kap in Kalist:
        for kdp in Kblist:
            Llist = [1,100,10000]
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
                    V =  Parameter(name='V', expression=1.66 * 10 ** (-13))
                    ka = Parameter(name = 'ka', expression = kap)
                    kd = Parameter(name='kd', expression=kdp)
                    k1 = Parameter(name='k1', expression=0.01)
                    k2 = Parameter(name='k2', expression=0.011)
                    d1 = Parameter(name='d1', expression=0.005)
                    d2 = Parameter(name='d2', expression=0.00022)
                    self.add_parameter([ka,kd,k1,k2,d1,d2,V])

                    # Define molecular species.
                    DNAoff = Species(name = 'DNAoff', initial_value = 1)
                    DNAon = Species(name = 'DNAon', initial_value = 0)
                    mRNA = Species(name='mRNA', initial_value=0)
                    Protein = Species(name='Protein', initial_value=0)
                    L = Species(name='L', initial_value=Lv)
                    self.add_species([DNAoff,DNAon,mRNA,Protein,L])

                    # Define reactions.
                    DNAoffandL_DNAon = Reaction(name = "DNAoffandL_DNAon", reactants = {DNAoff:1,L:1}, products = {DNAon:1},
                                  propensity_function = "ka*DNAoff*L*V")
                    DNAon_DNAoffandL = Reaction(name="DNAon_DNAoffandL", reactants={DNAon:1}, products={DNAoff:1,L:1},
                                        propensity_function="kd*DNAon")
                    DNAon_DNAonandmRNA = Reaction(name="DNAon_DNAonandmRNA", reactants={DNAon: 1}, products={DNAon: 1,mRNA:1},
                                        propensity_function="k1*DNAon")
                    mRNA_O = Reaction(name="mRNA_O", reactants={mRNA: 1}, products={},
                                        propensity_function="mRNA*d1")
                    mRNA_mRNAandProtein = Reaction(name="mRNA_mRNAandProtein", reactants={mRNA: 1},products={mRNA: 1, Protein: 1},
                                                  propensity_function="k2*mRNA")
                    Protein_O = Reaction(name="Protein_O", reactants={Protein: 1}, products={},
                                      propensity_function="Protein*d2")
                    self.add_reaction([DNAoffandL_DNAon,DNAon_DNAoffandL,DNAon_DNAonandmRNA,mRNA_O,mRNA_mRNAandProtein,Protein_O])
                    self.timespan(numpy.linspace(0, 10*1/0.00022,8000))
                mRNA_Synth_model = mRNA_Synth()
                results = []
                n_results = mRNA_Synth_model.run(solver = gillespy2.solvers.NumPySSASolver, number_of_trajectories=1,debug = False)
                print([kap,kdp,Lv])
                proteinfin = []
                for i in n_results:
                    proteinfin.append(i['Protein'][-1])
                input(proteinfin)
                proteinpoint = []
                for i in n_results:
                    proteinpoint.append(i['Protein'][200])
                print(n_results['time'][200])
                input(proteinpoint)
                proteinl = [0 for i in range(len(n_results[0]['time']))]
                mrna = []
                for i in n_results:
                    for j in range(len(i['Protein'])):
                        proteinl[j]+=i['Protein'][j]
                for i in range(len(proteinl)):
                    proteinl[i] /= len(n_results)
                plt.plot(n_results[0]['time']/3600,proteinl, label="Protein (Average)")
                plt.plot(n_results[0]['time']/3600,n_results[0]['Protein'],label = "Protein (Single Run)")
                plt.legend()
                plt.xlabel("time (hours)")
                plt.ylabel("count")
                plt.show()
                plt.plot(n_results[0]['time']/3600, n_results[0]['mRNA'],label = "mRNA")
                plt.plot(n_results[0]['time']/3600, n_results[0]['DNAon'], label="Activated DNA")
                plt.legend()
                plt.xlabel("time (hours)")
                plt.ylabel("count")
                plt.show()
                LRdist = []
                for i in n_results:
                    LRdist.append(i['Protein'][-1])
                    possibleLRlist.append(i['Protein'][-1])
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
            Infoperk[str([kap,kdp])] = information
            info[str([kap, kdp])] = [pLRgL, pLR, pLdict, pLRLdict]
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