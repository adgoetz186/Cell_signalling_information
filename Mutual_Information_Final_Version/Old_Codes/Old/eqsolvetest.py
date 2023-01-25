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

parameter_a_sweep = [-3.31,-3.06,-2.81,-2.56,-2.31]
parameter_b_sweep = [-3.5,-3,-2.5,-2,-1.5]
a = -2.81
b = -2.5
class EGFR_Signaling_Initial(Model):
    def __init__(self, parameter_values=None):
        # Initialize the model.
        Model.__init__(self, name="EGFR_Signaling_Initial")

        # Define parameters.
        k1b = Parameter(name='k1b', expression=10 ** (a))
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
        ckprod = Parameter(name='ckprod', expression=234423*10**b*10 ** (-3.6)/(10 ** (-3.6) + 10 ** (-3.56)))
        kndeg = Parameter(name='kndeg', expression=math.log(2)/((22.23)*(60**2)))
        knprod = Parameter(name='knprod', expression=21.58/(60**2))
        kmdeg = Parameter(name='kmdeg', expression=math.log(2)/((8.25)*(60**2)))
        kmprod = Parameter(name='kmprod', expression=math.log(2)/((8.25)*(60**2))*10.44)
        self.add_parameter([k1b, krec, ck2, kn1, ki, kn2, kp, kdp, ksi, ksrec, kdeg, ksdeg, ckprod,kndeg,knprod,kmdeg,kmprod])
        print(ckprod.expression)
        print(ki.expression)
        print(k1b.expression)
        initialVectorSolver = GEKKO()
        R,LR,LR2,sLR2,Ri,LRi,LR2i,sLR2i = [initialVectorSolver.Var(1) for i in range(8)]
        initialVectorSolver.Equations([ckprod.expression+Ri*krec.expression+LR*kn1.expression-R*ki.expression - R*k1b.expression==0,
                                       R*k1b.expression+LRi*krec.expression+LR2*kn2.expression-LR*kn1.expression-LR*ki.expression-LR*ck2.expression == 0,
                                       LR*ck2.expression+LR2i*krec.expression+sLR2*kdp.expression-LR2*kn2.expression-LR2*ki.expression-LR2*kp.expression== 0,
                                       LR2*kp.expression+sLR2i*ksrec.expression-sLR2*kdp.expression-sLR2*ksi.expression==0,
                                       R*ki.expression-Ri*kdeg.expression-Ri*krec.expression == 0,
                                       LR*ki.expression-LRi*kdeg.expression-LRi*krec.expression==0,
                                       LR2*ki.expression+sLR2i*kdp.expression-LR2i*krec.expression-LR2i*kp.expression-LR2i*kdeg.expression == 0,
                                       LR2i*kp.expression + sLR2*ksi.expression-sLR2i*ksdeg.expression-sLR2i*kdp.expression-sLR2i*ksrec.expression==0])
        initialVectorSolver.solve(disp=False)

        print([round(R.value[0]),isinstance(R.value[0],int),isinstance(round(R.value[0]),int)])

EGFR_model_Initial = EGFR_Signaling_Initial()