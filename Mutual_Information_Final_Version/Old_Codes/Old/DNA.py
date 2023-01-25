import numpy as np
np.max_line_width=np.inf
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
from scipy.stats import gamma
from scipy.stats import percentileofscore
from scipy.stats import kstest
from scipy.stats import nbinom


def homebrew_KS(observed,cdf_distribution,x_values):
    observed_cdf = np.array([percentileofscore(observed,x,kind='weak')/100 for x in x_values])
    plt.plot(x_values,observed_cdf)
    plt.plot(x_values,cdf_distribution)
    D = np.max(np.abs(observed_cdf - cdf_distribution))
    try1 = kstest(observed_cdf, cdf_distribution,N=100000)
    print(try1)
    try2 = kstest(observed_cdf, cdf_distribution,N=1000)
    print(try1)
    print(try2)
    print(D)
    plt.show()

def trial_receptor(x):
    Infoperk = {}
    info = {}
    #prove or discuss the ratios is all we need
    kstep = 0.5
    Kalist = [10**(3+i*kstep) for i in range(16)]
    Kblist = [10**(-5+i*kstep) for i in range(10)]
    Kalist = [10**11]
    Kblist = [10**(-3)]
    totalrunsneeded= len(Kblist)*len(Kalist)
    runcount = 0
    timestart = time.time()
    for kap in Kalist:
        for kdp in Kblist:
            Llist = [1,10,100,1000,10000]
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
                    self.timespan(np.linspace(0, 10*1/0.00022,80))
                mRNA_Synth_model = mRNA_Synth()

                class mRNA_Synth_nl(Model):
                  def __init__(self, parameter_values = None):
                    # Initialize the model.
                    Model.__init__(self, name = "mRNA_Synth_nl")

                    # Define parameters.
                    V =  Parameter(name='V', expression=1.66 * 10 ** (-13))
                    L = Parameter(name='L', expression=Lv)
                    ka = Parameter(name = 'ka', expression = kap)
                    kd = Parameter(name='kd', expression=kdp)
                    k1 = Parameter(name='k1', expression=0.01)
                    k2 = Parameter(name='k2', expression=0.011)
                    d1 = Parameter(name='d1', expression=0.005)
                    d2 = Parameter(name='d2', expression=0.00022)
                    self.add_parameter([ka,kd,k1,k2,d1,d2,V,L])

                    # Define molecular species.
                    DNAoff = Species(name = 'DNAoff', initial_value = 1)
                    DNAon = Species(name = 'DNAon', initial_value = 0)
                    mRNA = Species(name='mRNA', initial_value=0)
                    Protein = Species(name='Protein', initial_value=0)
                    self.add_species([DNAoff,DNAon,mRNA,Protein])

                    # Define reactions.
                    DNAoffandL_DNAon = Reaction(name = "DNAoffandL_DNAon", reactants = {DNAoff:1}, products = {DNAon:1},
                                  propensity_function = "ka*DNAoff*L*V")
                    DNAon_DNAoffandL = Reaction(name="DNAon_DNAoffandL", reactants={DNAon:1}, products={DNAoff:1},
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
                    self.timespan(np.linspace(0, 10*1/0.00022,80))
                mRNA_Synth_model_nl = mRNA_Synth_nl()
                results = []
                time_start = time.time()
                n_results_nl = mRNA_Synth_model_nl.run(solver=gillespy2.solvers.NumPySSASolver, number_of_trajectories=1000,
                                                 debug=False)
                n_results = mRNA_Synth_model.run(solver = gillespy2.solvers.NumPySSASolver, number_of_trajectories=1000,debug = False)
                print(time.time()-time_start)
                protein_at_final_time = []
                for i in n_results:
                    protein_at_final_time.append(i['Protein'][-1])
                protein_at_final_time = np.array(protein_at_final_time)
                continuous_range_of_values = np.linspace(np.min(protein_at_final_time)*.75,np.max(protein_at_final_time)*1.5,1000)
                discrete_range_of_values = np.arange(int(round(np.min(protein_at_final_time)*.75)),int(round(np.max(protein_at_final_time)*1.5)))
                tf = Lv
                k1 = kap
                k2 = kdp
                scale = 0.000120481927710843*(83.0*k1*tf + 500000000000000.0*k2)*(-2.7556e-22*k1**2*tf**2/(1.66e-13*k1*tf + k2)**2 + (1538735995700.0*k1**3*tf**3 + 9.5488429e+24*k1**2*k2*tf**2 + 4.8386758419e+22*k1**2*tf**2 + 1.682825e+36*k1*k2**2*tf + 2.01968465e+34*k1*k2*tf + 6.14243575e+31*k1*tf)/(149236407.0*k1**3*tf**3 + 2.6970435e+21*k1**2*k2*tf**2 + 4.69285569e+18*k1**2*tf**2 + 1.624725e+34*k1*k2**2*tf + 5.654043e+31*k1*k2*tf + 5.957325e+27*k1*tf + 3.2625e+46*k2**3 + 1.703025e+44*k2**2 + 3.58875e+40*k2))/(k1*tf)
                shape = 68890000.0*k1**2*tf**2*(4.0e-30/(1.66e-13*k1*tf + k2)**2)/(-2.7556e-22*k1**2*tf**2/(1.66e-13*k1*tf + k2)**2 + (1538735995700.0*k1**3*tf**3 + 9.5488429e+24*k1**2*k2*tf**2 + 4.8386758419e+22*k1**2*tf**2 + 1.682825e+36*k1*k2**2*tf + 2.01968465e+34*k1*k2*tf + 6.14243575e+31*k1*tf)/(149236407.0*k1**3*tf**3 + 2.6970435e+21*k1**2*k2*tf**2 + 4.69285569e+18*k1**2*tf**2 + 1.624725e+34*k1*k2**2*tf + 5.654043e+31*k1*k2*tf + 5.957325e+27*k1*tf + 3.2625e+46*k2**3 + 1.703025e+44*k2**2 + 3.58875e+40*k2))
                #distribution = gamma.pdf(range_of_values,31.5062379410726,scale = 3.11464113827816)
                distribution2 = gamma.pdf(continuous_range_of_values, shape, scale=scale)
                variance = shape*scale**2
                mean = shape*scale
                print(mean,variance)
                p = mean / variance
                n = mean * p / (1 - p)
                distribution3 = nbinom.pmf(discrete_range_of_values, n, p)
                #homebrew_KS(protein_at_final_time, gamma.cdf(range_of_values, shape, scale=scale), range_of_values)
                #plt.plot(range_of_values,distribution)
                plt.plot(discrete_range_of_values, distribution3)
                plt.plot(continuous_range_of_values,distribution2, '--')
                plt.hist(protein_at_final_time,bins=50,density=True)
                plt.title(f"Conditinoal probability p(R|S)\n s = {Lv}, or {round(Lv/(1.66 * 10 ** (-13))/(6.022*10**(23))*10**9,2)} nmol")
                plt.xlabel("Individual Protein Count")
                plt.ylabel("Probability")
                plt.legend(["Moment Prediction (Negative Binomial)", "Moment Prediction (Gamma)","Gillespie Simulation"])
                #plt.show()
                #plt.savefig(f"Final_Pres {Lv} g and nb")
                plt.clf()
                print(np.average(protein_at_final_time))
                print(np.var(protein_at_final_time))
                #input(n_results)
                print([kap,kdp,Lv])
                proteinl = [0 for i in range(len(n_results[0]['time']))]
                proteinlnl = [0 for i in range(len(n_results_nl[0]['time']))]
                mrna = []
                for i in n_results:
                    for j in range(len(i['Protein'])):
                        proteinl[j]+=i['Protein'][j]

                for i in n_results_nl:
                    for j in range(len(i['Protein'])):
                        proteinlnl[j]+=i['Protein'][j]

                for i in range(len(proteinl)):
                    proteinl[i] /= len(n_results)
                for i in range(len(proteinl)):
                    proteinlnl[i] /= len(n_results)
                plt.plot(n_results[0]['time'] / 3600, proteinl, label="True 2nd order reaction network")
                plt.plot(n_results[0]['time'] / 3600, proteinlnl, label="Network with 1st order assumption")
                plt.title("Comparison of 1st order assumption with true 2nd order rate\n1 DNA, high binding affinity")
                plt.xlabel("Time (Hours)")
                plt.ylabel("Protein Count")
                plt.legend()
                plt.show()
                #input(proteinl[-1])
                for i in range(len(n_results)):
                    plt.plot(n_results[0]['time'] / 3600, n_results[i]['Protein'], label="Protein (Average)")
                #plt.plot(n_results[0]['time']/3600,n_results[0]['Protein'],label = "Protein (Single Run)")
                #plt.legend()
                plt.title("10 Individual Trajectories For Single Ligand System")
                plt.xlabel("time (hours)")
                plt.ylabel("count")
                plt.show()
                #plt.close()
                plt.plot(n_results[0]['time']/3600, n_results[0]['mRNA'],label = "mRNA")
                plt.plot(n_results[0]['time']/3600, n_results[0]['DNAon'], label="Activated DNA")
                plt.legend()
                plt.xlabel("time (hours)")
                plt.ylabel("count")
                plt.show()
                #plt.close()
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