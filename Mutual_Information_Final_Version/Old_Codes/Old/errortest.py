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
EGFR_Initial_Results = {'time': [   0.,  600., 1200., 1800., 2400., 3000., 3600.], 'LR': [421., 406., 396., 387., 369., 344., 344.], 'LR2': [167., 178., 156., 131., 108., 116., 104.], 'LR2i': [182., 175., 167., 147., 123., 155., 133.], 'LRi': [139., 131., 139., 129., 128., 121., 108.], 'Nmkr': [8491., 8486., 8465., 8462., 8465., 8454., 8447.], 'R': [88002., 83597., 80235., 77777., 75818., 74422., 73156.], 'Ri': [29047., 28818., 28259., 27518., 26963., 26490., 25663.], 'mRNA': [12., 12., 12., 12., 12., 12., 12.], 'sLR2': [2509., 2278., 2168., 1991., 1879., 1763., 1810.], 'sLR2i': [2893., 2819., 2581., 2316., 2200., 2076., 2071.]}
a = -2.81
b = -2.5
fold = 2
class EGFR_Signaling_True(Model):
    def __init__(self, parameter_values=None):
        # Initialize the model.
        Model.__init__(self, name="EGFR_Signaling_True")

        # Define parameters.
        k1b = Parameter(name='k1b', expression=10 ** (a) * fold)
        krec = Parameter(name='krec', expression=10 ** (-3.56))
        ck2 = Parameter(name='ck2', expression=10 ** (-6.3))
        kn1 = Parameter(name='kn1', expression=10 ** (-0.52))
        ki = Parameter(name='ki', expression=10 ** (b) / 7223)
        kn2 = Parameter(name='kn2', expression=10 ** (-1.25))
        kp = Parameter(name='kp', expression=10 ** (0.05))
        kdp = Parameter(name='kdp', expression=10 ** (-1.15))
        ksi = Parameter(name='ksi', expression=10 ** (-2.5) / 7223)
        ksrec = Parameter(name='ksrec', expression=10 ** (-4.42))
        kdeg = Parameter(name='kdeg', expression=10 ** (-3.6))
        ksdeg = Parameter(name='ksdeg', expression=10 ** (-2.5))
        ckprod = Parameter(name='ckprod', expression=16)
        print(ckprod.expression)
        kndeg = Parameter(name='kndeg', expression=math.log(2) / ((22.23) * (60 ** 2)))
        knprod = Parameter(name='knprod', expression=21.58 / (60 ** 2))
        kmdeg = Parameter(name='kmdeg', expression=math.log(2) / ((8.25) * (60 ** 2)))
        kmprod = Parameter(name='kmprod', expression=math.log(2) / ((8.25) * (60 ** 2)) * 10.44)
        self.add_parameter(
            [k1b, krec, ck2, kn1, ki, kn2, kp, kdp, ksi, ksrec, kdeg, ksdeg, ckprod, kndeg, knprod, kmdeg, kmprod])

        # Define molecular species.
        R = Species(name='R', initial_value=EGFR_Initial_Results['R'][-1])
        LR = Species(name='LR', initial_value=EGFR_Initial_Results['LR'][-1])
        LR2 = Species(name='LR2', initial_value=EGFR_Initial_Results['LR2'][-1])
        sLR2 = Species(name='sLR2', initial_value=EGFR_Initial_Results['sLR2'][-1])
        Ri = Species(name='Ri', initial_value=EGFR_Initial_Results['Ri'][-1])
        LRi = Species(name='LRi', initial_value=EGFR_Initial_Results['LRi'][-1])
        LR2i = Species(name='LR2i', initial_value=EGFR_Initial_Results['LR2i'][-1])
        sLR2i = Species(name='sLR2i', initial_value=EGFR_Initial_Results['sLR2i'][-1])
        Nmkr = Species(name='Nmkr', initial_value=EGFR_Initial_Results['Nmkr'][-1])
        mRNA = Species(name='mRNA', initial_value=EGFR_Initial_Results['mRNA'][-1])
        self.add_species([R, LR, LR2, sLR2, Ri, LRi, LR2i, sLR2i, Nmkr, mRNA])

        # Define reactions.
        R_LR = Reaction(name="R_LR", reactants={R: 1}, products={LR: 1},
                        propensity_function="k1b*R")
        LR_R = Reaction(name="LR_R", reactants={LR: 1}, products={R: 1},
                        propensity_function="kn1*LR")
        O_R = Reaction(name="O_R", reactants={}, products={R: 1},
                       propensity_function="ckprod")
        RandLR_LR2 = Reaction(name="RandLR_LR2", reactants={R: 1, LR: 1}, products={LR2: 1},
                              propensity_function="ck2*R*LR")
        LR2_RandLR = Reaction(name="LR2_RandLR", reactants={LR2: 1}, products={R: 1, LR: 1},
                              propensity_function="kn2*LR2")
        LR2_sLR2 = Reaction(name="LR2_sLR2", reactants={LR2: 1}, products={sLR2: 1},
                            propensity_function="kp*LR2")
        sLR2_LR2 = Reaction(name="sLR2_LR2", reactants={sLR2: 1}, products={LR2: 1},
                            propensity_function="kdp*sLR2")
        LR2i_sLR2i = Reaction(name="LR2i_sLR2i", reactants={LR2i: 1}, products={sLR2i: 1},
                              propensity_function="kp*LR2i")
        sLR2i_LR2i = Reaction(name="sLR2i_LR2i", reactants={sLR2i: 1}, products={LR2i: 1},
                              propensity_function="kdp*sLR2i")
        RandNmkr_RiandNmkr = Reaction(name="RandNmkr_RiandNmkr", reactants={R: 1, Nmkr: 1}, products={Ri: 1, Nmkr: 1},
                                      propensity_function="ki*R*Nmkr")
        LRandNmkr_LRiandNmkr = Reaction(name="LRandNmkr_LRiandNmkr", reactants={LR: 1, Nmkr: 1},
                                        products={LRi: 1, Nmkr: 1},
                                        propensity_function="ki*LR*Nmkr")
        LR2andNmkr_LR2iandNmkr = Reaction(name="LR2andNmkr_LR2iandNmkr", reactants={LR2: 1, Nmkr: 1},
                                          products={LR2i: 1, Nmkr: 1},
                                          propensity_function="ki*LR2*Nmkr")
        sLR2andNmkr_sLR2iandNmkr = Reaction(name="sLR2andNmkr_sLR2iandNmkr", reactants={sLR2: 1, Nmkr: 1},
                                            products={sLR2i: 1, Nmkr: 1},
                                            propensity_function="ksi*sLR2*Nmkr")
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
        mRNA_mRNAandNmkr = Reaction(name="mRNA_mRNAandNmkr", reactants={mRNA: 1}, products={Nmkr: 1, mRNA: 1},
                                    propensity_function="knprod*mRNA")
        Nmkr_O = Reaction(name="Nmkr_O", reactants={Nmkr: 1}, products={},
                          propensity_function="kndeg*Nmkr")
        mRNA_O = Reaction(name="mRNA_O", reactants={mRNA: 1}, products={},
                          propensity_function="kmdeg*mRNA")
        O_mRNA = Reaction(name="O_mRNA", reactants={}, products={mRNA: 1},
                          propensity_function="kmprod")
        self.add_reaction(
            [R_LR, LR_R, O_R, RandLR_LR2, LR2_RandLR, LR2_sLR2, sLR2_LR2, LR2i_sLR2i, sLR2i_LR2i, RandNmkr_RiandNmkr,
             LRandNmkr_LRiandNmkr, LR2andNmkr_LR2iandNmkr, sLR2andNmkr_sLR2iandNmkr, Ri_R, LRi_LR, LR2i_LR2, sLR2i_sLR2,
             Ri_O, LRi_O,
             LR2i_O, sLR2i_O, mRNA_mRNAandNmkr, Nmkr_O, mRNA_O, O_mRNA])
        self.timespan(numpy.linspace(0, 1800, 1801))


EGFR_model_True = EGFR_Signaling_True()
#EGFR_True_Results = EGFR_model_True.run(solver=gillespy2.solvers.NumPySSASolver, number_of_trajectories=1, seed=None,debug=False, show_labels=True)[0]
EGFR_True_Results = EGFR_model_True.run(solver=gillespy2.solvers.ODESolver, number_of_trajectories=1, seed=None,debug=False, show_labels=True)[0]
print(EGFR_True_Results)
plt.plot(EGFR_True_Results["time"],EGFR_True_Results['sLR2'])
plt.show()