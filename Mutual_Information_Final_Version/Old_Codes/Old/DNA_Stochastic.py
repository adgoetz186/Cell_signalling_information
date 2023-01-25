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
import numpy as np


class mRNA_Synth(Model):
    def __init__(self, parameter_values=None):
        Model.__init__(self, name="A_B")
        k1 = Parameter(name='k1', expression=5)
        k2 = Parameter(name='k2', expression=1)
        self.add_parameter([k1,k2])

        # Define molecular species.
        A = Species(name = 'A', initial_value = 1000000)
        B = Species(name = 'B', initial_value = 0)
        self.add_species([A,B])

        # Define reactions.
        A_B = Reaction(name = "A_B", reactants = {A:1}, products = {B:1},
                      propensity_function = "k1*A")
        B_A = Reaction(name="B_A", reactants={B:1}, products={A:1},
                            propensity_function="k2*B")
        self.add_reaction([A_B,B_A])
        self.timespan(numpy.linspace(0, 3,100))
mRNA_Synth_model = mRNA_Synth()
results = []
n_results = mRNA_Synth_model.run(solver = gillespy2.solvers.NumPySSASolver, number_of_trajectories=1,debug = False)
n_results_det = mRNA_Synth_model.run(solver = gillespy2.solvers.ODESolver)
total_A_array = np.zeros((np.shape(n_results[0]['A'])))
total_B_array = np.zeros((np.shape(n_results[0]['B'])))
for i in n_results:
    total_A_array+= n_results['A']
    total_B_array += n_results['B']
total_A_array = total_A_array/len(n_results)
total_B_array = total_B_array/len(n_results)
plt.plot(n_results_det[0]['time'], n_results_det[0]['A'])
plt.plot(n_results_det[0]['time'], n_results_det[0]['B'])
plt.plot(n_results[0]['time'], total_A_array,'--')
plt.plot(n_results[0]['time'], total_B_array,'--')
plt.title("Gillespie Algorithm Compared to Deterministic Model")
plt.xlabel("Time (s)")
plt.ylabel("Molecule Count")
plt.legend(["Molecule A (Deterministic)","Molecule B (Deterministic)","Molecule A (Gillespie)","Molecule B (Gillespie)"])
plt.show()
print(n_results_det)