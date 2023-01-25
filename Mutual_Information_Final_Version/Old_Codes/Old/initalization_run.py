import numpy
import gillespy2
from gillespy2 import Model, Species, Parameter, Reaction
from gillespy2.solvers.numpy.basic_ode_solver import BasicODESolver
import matplotlib.pyplot as plt
import time
import math
import os
from multiprocessing import Pool
import multiprocessing
class EGFR_Signaling(Model):
  def __init__(self, parameter_values = None):
    # Initialize the model.
    Model.__init__(self, name = "EGFR_Signaling")

    # Define parameters.
    k1b = Parameter(name = 'k1b', expression = 10**(-2.81))
    krec = Parameter(name='krec', expression=10 ** (-3.56))
    ck2 = Parameter(name='ck2', expression=10 ** (-6.3))
    kn1 = Parameter(name = 'kn1', expression = 10 ** (-0.52))
    ki = Parameter(name='ki', expression=10 ** (-3.83)/7223)
    kn2 = Parameter(name='kn2', expression=10 ** (-1.25))
    kp = Parameter(name='kp', expression=10 ** (0.05))
    kdp = Parameter(name = 'kdp',   expression = 10 ** (-1.15))
    ksi = Parameter(name='ksi', expression=10 ** (-2.5)/7223)
    ksrec = Parameter(name='ksrec', expression=10 ** (-4.42))
    kdeg = Parameter(name='kdeg', expression=10 ** (-3.6))
    ksdeg = Parameter(name='ksdeg', expression=10 ** (-2.5))
    ckprod = Parameter(name='ckprod', expression=16.54)
    self.add_parameter([k1b, krec,ck2,kn1,ki,kn2,kp,kdp,ksi,ksrec,kdeg,ksdeg,ckprod])

    # Define molecular species.
    R = Species(name = 'R', initial_value = 234423)
    LR = Species(name = 'LR', initial_value = 0)
    LR2 = Species(name='LR2', initial_value=0)
    sLR2 = Species(name='sLR2', initial_value=0)
    Ri = Species(name='Ri', initial_value=65843)
    LRi = Species(name='LRi', initial_value=0)
    LR2i = Species(name='LR2i', initial_value=0)
    sLR2i = Species(name='sLR2i', initial_value=0)
    C = Species(name='C', initial_value=7223)
    self.add_species([R,LR,LR2,sLR2,Ri,LRi,LR2i,sLR2i,C])

    # Define reactions.
    R_LR = Reaction(name = "R_LR", reactants = {R:1}, products = {LR:1},
                  propensity_function = "k1b*R")
    LR_R = Reaction(name="LR_R", reactants={LR: 1}, products={R: 1},
                    propensity_function="kn1*LR")
    O_R = Reaction(name = "O_R", reactants = {}, products = {R:1},
                  propensity_function = "ckprod")
    RandLR_LR2 = Reaction(name="RandLR_LR2", reactants={R:1,LR:1}, products={LR2: 1},
                   propensity_function="ck2*R*LR")
    LR2_RandLR = Reaction(name="LR2_RandLR", reactants={LR2: 1}, products={R:1,LR:1},
                          propensity_function="kn2*LR2")
    LR2_sLR2 = Reaction(name="LR2_sLR2", reactants={LR2: 1}, products={sLR2: 1},
                  propensity_function="kp*LR2")
    sLR2_LR2 = Reaction(name="sLR2_LR2", reactants={sLR2: 1}, products={LR2: 1},
                  propensity_function="kdp*sLR2")
    LR2i_sLR2i = Reaction(name="LR2i_sLR2i", reactants={LR2i: 1}, products={sLR2i: 1},
                        propensity_function="kp*LR2i")
    sLR2i_LR2i = Reaction(name="sLR2i_LR2i", reactants={sLR2i: 1}, products={LR2i: 1},
                          propensity_function="kdp*sLR2i")
    RandC_RiandC = Reaction(name="RandC_RiandC", reactants={R: 1,C:1}, products={Ri: 1,C:1},
                          propensity_function="ki*R*C")
    LRandC_LRiandC = Reaction(name="LRandC_LRiandC", reactants={LR: 1, C: 1}, products={LRi: 1, C: 1},
                            propensity_function="ki*LR*C")
    LR2andC_LR2iandC = Reaction(name="LR2andC_LR2iandC", reactants={LR2: 1, C: 1}, products={LR2i: 1, C: 1},
                              propensity_function="ki*LR2*C")
    sLR2andC_sLR2iandC = Reaction(name="sLR2andC_sLR2iandC", reactants={sLR2: 1, C: 1}, products={sLR2i: 1, C: 1},
                                propensity_function="ksi*sLR2*C")
    Ri_R = Reaction(name="Ri_R", reactants={Ri: 1}, products={R: 1},
                                propensity_function="krec*Ri")
    LRi_LR = Reaction(name="LRi_LR", reactants={LRi: 1}, products={LR: 1},
                    propensity_function="krec*LRi")
    LR2i_LR2 = Reaction(name="LR2i_LR2", reactants={LR2i: 1}, products={LR2: 1},
                    propensity_function="krec*LR2i")
    sLR2i_sLR2 = Reaction(name="sLR2i_sLR2", reactants={sLR2i: 1}, products={sLR2: 1},
                    propensity_function="ksrec*sLR2i")
    Ri_O = Reaction(name="Ri_O", reactants={Ri: 1}, products={},
                          propensity_function="kdeg*Ri")
    LRi_O = Reaction(name="LRi_O", reactants={LRi: 1}, products={},
                    propensity_function="kdeg*LRi")
    LR2i_O = Reaction(name="LR2i_O", reactants={LR2i: 1}, products={},
                     propensity_function="kdeg*LR2i")
    sLR2i_O = Reaction(name="sLR2i_O", reactants={sLR2i: 1}, products={},
                      propensity_function="ksdeg*sLR2i")
    self.add_reaction([R_LR,LR_R,O_R,RandLR_LR2,LR2_RandLR,LR2_sLR2,sLR2_LR2,LR2i_sLR2i,sLR2i_LR2i,RandC_RiandC,LRandC_LRiandC,LR2andC_LR2iandC,sLR2andC_sLR2iandC,Ri_R,LRi_LR,LR2i_LR2,sLR2i_sLR2,Ri_O,LRi_O,LR2i_O,sLR2i_O])
    self.timespan(numpy.linspace(0, 100, 101))
EGFR_model = EGFR_Signaling()
for i in range(0):
    print(i)
    class mRNA_Synth(Model):
      def __init__(self, parameter_values = None):
        # Initialize the model.
        Model.__init__(self, name = "mass_action")

        # Define parameters.
        kmprod = Parameter(name = 'kmprod', expression = math.log(2)/((8.25)*(60**2))*10.44)
        kmdeg = Parameter(name='kmdeg', expression=math.log(2)/((8.25)*(60**2)))
        knprod = Parameter(name='knprod', expression= 21.58/(60**2))
        kndeg = Parameter(name = 'kndeg', expression = math.log(2)/((22.23)*(60**2)))
        self.add_parameter([kmprod,kmdeg,knprod,kndeg])

        # Define molecular species.
        mRNA = Species(name = 'mRNA', initial_value = 0)
        Nmkr = Species(name = 'Nmkr', initial_value = 0)
        self.add_species([mRNA,Nmkr])

        # Define reactions.
        O_mRNA = Reaction(name = "O_mRNA", reactants = {}, products = {mRNA:1},
                      propensity_function = "kmprod")
        mRNA_O = Reaction(name="mRNA_O", reactants={mRNA: 1}, products={},
                        propensity_function="kmdeg*mRNA")
        mRNA_mRNAandNmkr = Reaction(name = "mRNA_mRNAandNmkr", reactants = {mRNA:1}, products = {mRNA:1,Nmkr:1},
                      propensity_function = "knprod*mRNA")
        Nmkr_O = Reaction(name="Nmkr_O", reactants={Nmkr:1}, products={},
                       propensity_function="kndeg*Nmkr")
        self.add_reaction([O_mRNA,mRNA_O,mRNA_mRNAandNmkr,Nmkr_O])
        #self.timespan(numpy.linspace(0, 100, 101))
    mRNA_Synth_model = mRNA_Synth()

    timeStart = time.time()
    d_results = EGFR_model.run(solver = gillespy2.solvers.TauLeapingSolver, t=10, number_of_trajectories=1, seed=None, debug=False, show_labels=True)
    print("Time Tau = " + str(time.time() - timeStart))
    timeStart = time.time()
    n_results = mRNA_Synth_model.run(solver = gillespy2.solvers.NumPySSASolver, t=1728000, number_of_trajectories=1, seed=None, debug=False, show_labels=True)
    print("Time NM = " + str(time.time() - timeStart))
    plt.plot(n_results[0]['time'], n_results[0]['Nmkr'], '-r', label='U')
    plt.plot(n_results[0]['time'], n_results[0]['mRNA'], '-r', label='U')
    plt.title('Stochastic Switch')
    plt.legend(loc = 'best')