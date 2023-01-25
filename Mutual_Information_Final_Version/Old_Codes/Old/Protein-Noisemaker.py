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


def Noise_Maker(x):


    for fold in range(2,4):
        parameter_a_sweep = [-3.31,-2.81,-2.31]
        parameter_b_sweep = [-3,-2.5,-2]
        parameter_a_sweep = [-2.81]
        parameter_b_sweep = [-2.5]
        resultsHolder = []

        for a in parameter_a_sweep:
            for b in parameter_b_sweep:
                resultsHoldersublist = [str(a),str(b)]
                timeResults = []
                timebegin = time.time()
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
                    mRNA = Species(name = 'mRNA', initial_value = 10)
                    Nmkr = Species(name = 'Nmkr', initial_value = 7225)
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
                    self.timespan(numpy.linspace(0, 2000000,6))
                kmprod = math.log(2) / ((8.25) * (60 ** 2)) * 10.44
                kmdeg = math.log(2) / ((8.25) * (60 ** 2))
                knprod = 21.58 / (60 ** 2)
                kndeg = math.log(2) / ((22.23) * (60 ** 2))
                ckproddd = 16
                Nmkrss = kmprod*knprod/(kmdeg*kndeg)
                mRNA_Synth_model = mRNA_Synth()
                results = []
                n_results = mRNA_Synth_model.run(solver = gillespy2.solvers.NumPySSASolver, number_of_trajectories=1, seed=None, debug=False, show_labels=True)
                results.append([n_results[0]['mRNA'][-1],n_results[0]['Nmkr'][-1]])
                timeResults.append(time.time()-timebegin)
                print("Time first: " + str(timeResults[-1]))
                timebegin = time.time()
                # Below we assume ckprod is the same for each cell
                # If we wanted to have it determined for each cell, divide ki by the average noisemaker then
                # multiply ki in the solve equations by the actual noisemaker count of the system
                NoLigandSolver = GEKKO(remote=False)
                krec = 10 ** (-3.56)
                kdeg = 10 ** (-3.6)
                ki = 10 ** (-3.83)
                R, Ri, ckprod = [NoLigandSolver.Var(1) for i in range(3)]
                NoLigandSolver.Equations([R + Ri == 3 * 10 ** 5,
                                          ckprod + Ri * krec - R * ki == 0,
                                          R * ki - Ri * krec - Ri * kdeg == 0])
                NoLigandSolver.solve(disp=False)
                print("ckprod = " + str(ckprod.value[0]))
                print("R = " + str(R.value[0]))
                print("Ri = " + str(Ri.value[0]))
                ckprodd = ckprod.value[0]
                print(ckprodd)
                class EGFR_Signaling_Initial(Model):
                    def __init__(self, parameter_values=None):
                        # Initialize the model.
                        Model.__init__(self, name="EGFR_Signaling_Initial")




                        # Define parameters.
                        k1b = Parameter(name='k1b', expression=10 ** (a))
                        krec = Parameter(name='krec', expression=10 ** (-3.56))
                        ck2 = Parameter(name='ck2', expression=10 ** (-6.3))
                        kn1 = Parameter(name='kn1', expression=10 ** (-0.52))
                        ki = Parameter(name='ki', expression=10 ** (-3.83)/Nmkrss)
                        kn2 = Parameter(name='kn2', expression=10 ** (-1.25))
                        kp = Parameter(name='kp', expression=10 ** (0.05))
                        kdp = Parameter(name='kdp', expression=10 ** (-1.15))
                        ksi = Parameter(name='ksi', expression=10 ** (b) /Nmkrss)
                        ksrec = Parameter(name='ksrec', expression=10 ** (-4.42))
                        kdeg = Parameter(name='kdeg', expression=10 ** (-3.6))
                        ksdeg = Parameter(name='ksdeg', expression=10 ** (-2.5))
                        ckprod = Parameter(name='ckprod', expression=ckprodd)
                        print(ckprod.expression)
                        print(n_results[0]['Nmkr'][-1])
                        kndeg = Parameter(name='kndeg', expression=math.log(2)/((22.23)*(60**2)))
                        knprod = Parameter(name='knprod', expression=21.58/(60**2))
                        kmdeg = Parameter(name='kmdeg', expression=math.log(2)/((8.25)*(60**2)))
                        kmprod = Parameter(name='kmprod', expression=math.log(2)/((8.25)*(60**2))*10.44)
                        self.add_parameter([k1b, krec, ck2, kn1, ki, kn2, kp, kdp, ksi, ksrec, kdeg, ksdeg, ckprod,kndeg,knprod,kmdeg,kmprod])

                        print("solve start")
                        solvetime = time.time()
                        initialVectorSolver = GEKKO(remote=False)
                        R,LR,LR2,sLR2,Ri,LRi,LR2i,sLR2i = [initialVectorSolver.Var(1) for i in range(8)]
                        initialVectorSolver.Equations([ckprod.expression+Ri*krec.expression+LR*kn1.expression-R*ki.expression*n_results[0]['Nmkr'][-1] - R*k1b.expression==0,
                                           R*k1b.expression+LRi*krec.expression+LR2*kn2.expression-LR*kn1.expression-LR*ki.expression*n_results[0]['Nmkr'][-1]-LR*ck2.expression*R == 0,
                                           LR*ck2.expression*R+LR2i*krec.expression+sLR2*kdp.expression-LR2*kn2.expression-LR2*ki.expression*n_results[0]['Nmkr'][-1]-LR2*kp.expression== 0,
                                           LR2*kp.expression+sLR2i*ksrec.expression-sLR2*kdp.expression-sLR2*ksi.expression*n_results[0]['Nmkr'][-1]==0,
                                           R*ki.expression*n_results[0]['Nmkr'][-1]-Ri*kdeg.expression-Ri*krec.expression == 0,
                                           LR*ki.expression*n_results[0]['Nmkr'][-1]-LRi*kdeg.expression-LRi*krec.expression==0,
                                           LR2*ki.expression*n_results[0]['Nmkr'][-1]+sLR2i*kdp.expression-LR2i*krec.expression-LR2i*kp.expression-LR2i*kdeg.expression == 0,
                                           LR2i*kp.expression + sLR2*ksi.expression*n_results[0]['Nmkr'][-1]-sLR2i*ksdeg.expression-sLR2i*kdp.expression-sLR2i*ksrec.expression==0])
                        initialVectorSolver.solve(disp=False)
                        print("solve time = " + str(time.time()-solvetime))
                        print([R.value[0],LR.value[0],LR2.value[0],sLR2.value[0],Ri.value[0],LRi.value[0],LR2i.value[0],sLR2i.value[0]])
                        # Define molecular species.
                        R = Species(name='R', initial_value=round(R.value[0]))
                        LR = Species(name='LR', initial_value=round(LR.value[0]))
                        LR2 = Species(name='LR2', initial_value=round(LR2.value[0]))
                        sLR2 = Species(name='sLR2', initial_value=round(sLR2.value[0]))
                        Ri = Species(name='Ri', initial_value=round(Ri.value[0]))
                        LRi = Species(name='LRi', initial_value=round(LRi.value[0]))
                        LR2i = Species(name='LR2i', initial_value=round(LR2i.value[0]))
                        sLR2i = Species(name='sLR2i', initial_value=round(sLR2i.value[0]))
                        Nmkr = Species(name='Nmkr', initial_value=n_results[0]['Nmkr'][-1])
                        mRNA = Species(name='mRNA', initial_value=n_results[0]['mRNA'][-1])
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
                        LRandNmkr_LRiandNmkr = Reaction(name="LRandNmkr_LRiandNmkr", reactants={LR: 1, Nmkr: 1}, products={LRi: 1, Nmkr: 1},
                                                  propensity_function="ki*LR*Nmkr")
                        LR2andNmkr_LR2iandNmkr = Reaction(name="LR2andNmkr_LR2iandNmkr", reactants={LR2: 1, Nmkr: 1}, products={LR2i: 1, Nmkr: 1},
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
                        mRNA_mRNAandNmkr = Reaction(name="mRNA_mRNAandNmkr", reactants={mRNA:1}, products={Nmkr:1,mRNA:1},
                                           propensity_function="knprod*mRNA")
                        Nmkr_O = Reaction(name="Nmkr_O", reactants={Nmkr: 1}, products={},
                                          propensity_function="kndeg*Nmkr")
                        mRNA_O = Reaction(name="mRNA_O", reactants={mRNA: 1}, products={},
                                                    propensity_function="kmdeg*mRNA")
                        O_mRNA = Reaction(name="O_mRNA", reactants={}, products={mRNA: 1},
                                          propensity_function="kmprod")
                        self.add_reaction(
                            [R_LR, LR_R, O_R, RandLR_LR2, LR2_RandLR, LR2_sLR2, sLR2_LR2, LR2i_sLR2i, sLR2i_LR2i, RandNmkr_RiandNmkr,
                             LRandNmkr_LRiandNmkr, LR2andNmkr_LR2iandNmkr, sLR2andNmkr_sLR2iandNmkr, Ri_R, LRi_LR, LR2i_LR2, sLR2i_sLR2, Ri_O, LRi_O,
                             LR2i_O, sLR2i_O,mRNA_mRNAandNmkr,Nmkr_O,mRNA_O,O_mRNA])
                        self.timespan(numpy.linspace(0, 3600, 7))
                EGFR_model_Initial = EGFR_Signaling_Initial()
                EGFR_Initial_Results = EGFR_model_Initial.run(solver=gillespy2.solvers.TauLeapingSolver, number_of_trajectories=1, seed=None,
                                                 debug=False, show_labels=True)[0]
                print(EGFR_Initial_Results)
                timeResults.append(time.time()-timebegin)
                print("Time Mid: " + str(timeResults[-1]))
                timebegin = time.time()
                foldResults = []
                class EGFR_Signaling_True(Model):
                    def __init__(self, parameter_values=None):
                        # Initialize the model.
                        Model.__init__(self, name="EGFR_Signaling_True")

                        # Define parameters.
                        k1b = Parameter(name='k1b', expression=10 ** (a)*fold)
                        krec = Parameter(name='krec', expression=10 ** (-3.56))
                        ck2 = Parameter(name='ck2', expression=10 ** (-6.3))
                        kn1 = Parameter(name='kn1', expression=10 ** (-0.52))
                        ki = Parameter(name='ki', expression=10 ** (-3.83)/Nmkrss)
                        kn2 = Parameter(name='kn2', expression=10 ** (-1.25))
                        kp = Parameter(name='kp', expression=10 ** (0.05))
                        kdp = Parameter(name='kdp', expression=10 ** (-1.15))
                        ksi = Parameter(name='ksi', expression=10 ** (b) /Nmkrss)
                        ksrec = Parameter(name='ksrec', expression=10 ** (-4.42))
                        kdeg = Parameter(name='kdeg', expression=10 ** (-3.6))
                        ksdeg = Parameter(name='ksdeg', expression=10 ** (-2.5))
                        ckprod = Parameter(name='ckprod', expression=ckprodd)
                        kndeg = Parameter(name='kndeg', expression=math.log(2)/((22.23)*(60**2)))
                        knprod = Parameter(name='knprod', expression=21.58/(60**2))
                        kmdeg = Parameter(name='kmdeg', expression=math.log(2)/((8.25)*(60**2)))
                        kmprod = Parameter(name='kmprod', expression=math.log(2)/((8.25)*(60**2))*10.44)
                        self.add_parameter([k1b, krec, ck2, kn1, ki, kn2, kp, kdp, ksi, ksrec, kdeg, ksdeg, ckprod,kndeg,knprod,kmdeg,kmprod])

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
                        LRandNmkr_LRiandNmkr = Reaction(name="LRandNmkr_LRiandNmkr", reactants={LR: 1, Nmkr: 1}, products={LRi: 1, Nmkr: 1},
                                                  propensity_function="ki*LR*Nmkr")
                        LR2andNmkr_LR2iandNmkr = Reaction(name="LR2andNmkr_LR2iandNmkr", reactants={LR2: 1, Nmkr: 1}, products={LR2i: 1, Nmkr: 1},
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
                        mRNA_mRNAandNmkr = Reaction(name="mRNA_mRNAandNmkr", reactants={mRNA:1}, products={Nmkr:1,mRNA:1},
                                           propensity_function="knprod*mRNA")
                        Nmkr_O = Reaction(name="Nmkr_O", reactants={Nmkr: 1}, products={},
                                          propensity_function="kndeg*Nmkr")
                        mRNA_O = Reaction(name="mRNA_O", reactants={mRNA: 1}, products={},
                                                    propensity_function="kmdeg*mRNA")
                        O_mRNA = Reaction(name="O_mRNA", reactants={}, products={mRNA: 1},
                                          propensity_function="kmprod")
                        self.add_reaction(
                            [R_LR, LR_R, O_R, RandLR_LR2, LR2_RandLR, LR2_sLR2, sLR2_LR2, LR2i_sLR2i, sLR2i_LR2i, RandNmkr_RiandNmkr,
                             LRandNmkr_LRiandNmkr, LR2andNmkr_LR2iandNmkr, sLR2andNmkr_sLR2iandNmkr, Ri_R, LRi_LR, LR2i_LR2, sLR2i_sLR2, Ri_O, LRi_O,
                             LR2i_O, sLR2i_O,mRNA_mRNAandNmkr,Nmkr_O,mRNA_O,O_mRNA])
                        self.timespan(numpy.linspace(0, 1800, 61))

                EGFR_model_True = EGFR_Signaling_True()
                EGFR_True_Results = EGFR_model_True.run(solver=gillespy2.solvers.NumPySSASolver, number_of_trajectories=1, seed=None,
                                                   debug=False, show_labels=True)[0]
                print(EGFR_True_Results)
                print("test")
                timeResults.append(time.time()-timebegin)
                print("Time Last: " + str(timeResults[-1]))
                foldResults.append(EGFR_True_Results)
                print(timeResults)
                print(foldResults)
                print(len(foldResults))
                print("testy")
                resultsHoldersublist.append(foldResults)
                resultsHolder.append(resultsHoldersublist)
        for i in resultsHolder:
            filename = "Fold_is_" + str(fold) +  "____a_is_" + i[0]+"____b_is_"+i[1]
            if os.path.exists(filename):
                with open(filename, "a+") as runDataFile:
                    runDataFile.write(str(i[2]).replace("\n","").replace("  ","") + "\n")
            else:
                with open(filename, "w+") as runDataFile:
                    runDataFile.write(str(i[2]).replace("\n","").replace("  ","") + "\n")
if __name__ == '__main__':
    #os.chdir('/blue/pdixit/agoetz/Raw_Output/')
    timestart = time.time()
    #Noise_Maker(1)
    for repeat in range(5):
        with Pool(multiprocessing.cpu_count()) as p:
            p.map(Noise_Maker,range(multiprocessing.cpu_count()))
            print(1)
    print(time.time()-timestart)