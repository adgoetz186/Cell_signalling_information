from Mutual_Information_Final_Version.functions.Models.IGFR_Akt_FoxO import model
import numpy as np
from pysb.simulator import BngSimulator
import multiprocessing as mp
import functools
import Mutual_Information_Final_Version.functions.file_processing_functions.user_input_functions as uif
import time
import random
import matplotlib.pyplot as plt
import copy as cp
plt.rcParams.update({'font.size': 15})

# This program generates outcomes of the IGFR/Akt/FoxO pathway using the Gillespie algorithm
# These results are used to test the moment closure accuracy





# _____ Load arguments BEGIN _____
# Specifies path to file which specifies user specified arguments
user_input_file_path = "/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Figure_Generating_Programs/S_Gillespie_Supplimentary/General_Data_Generation/IGFR/ssa_IGFR_args.txt"
user_arguments = uif.read_file_to_arg_dict(user_input_file_path)
# _____ Load arguments END _____

# _____ File path declarations BEGIN _____
# Specifies path to drawn model cell parameters
params_file_path = user_arguments["params_file_path"]
# Specifies path to output directory
output_folder = user_arguments["output_dir"]
# _____ File path declarations END _____

# _____ Loading files BEGIN _____
# Loads parameters of drawn cells
params = 10 ** np.loadtxt(params_file_path, delimiter=",")
# _____ Loading files END _____

# _____ Main code BEGIN _____
# time allowed before sampling cells with no dose
t_f1 = 3600
t_f1_array = np.linspace(0, t_f1, 2)

# time points of main simulation
t_main = np.arange(0, 180*61, 180)
# list of doses
L = [0,0.01,0.015,0.020,0.025,0.050]
dose_colors = ["black","indigo","blue","green","orange","red"]
# number of cells to use
cell_count = user_arguments["cell_count"]

# The set of parameters used elsewhere doesnt allow for non equilibrium initial conditions
# this transforms the parameters into a set which can allow for a relevant subset of
# non equilibrium initial conditions
kprod = params[:, 0] * params[:, 1]
akt_0 = params[:, 12].astype(int)
n_foxo_0 = (params[:, 13] * params[:, 10] / (params[:, 10] + params[:, 11])).astype(int)
c_foxo_0 = (params[:, 13] * params[:, 11] / (params[:, 10] + params[:, 11])).astype(int) + 1
igfr_0 = params[:, 0].astype(int)
new_params = np.zeros((np.shape(params)[0], np.shape(params)[1] + 3))
new_params[:, 0] = params[:, 0] * params[:, 1]
new_params[:, 1:12] = params[:, 1:12]
new_params[:, 12] = igfr_0
new_params[:, 13] = akt_0
new_params[:, 14] = c_foxo_0
new_params[:, 15] = n_foxo_0
time_start = time.time()

# Runs the ssa algorithm
for cell_ind in range(cell_count):
	cell_run_array = np.zeros((len(L),np.size(t_main)))
	total_cell = np.shape(new_params)[0]
	sim = BngSimulator(model, tspan=t_f1_array, param_values=new_params[cell_ind])
	# The seed is incredibly important for parallel processing since the time is used to create random seeds
	# often all parallel workers will end up with the same time which makes them all give identical results
	# this prevents workers from having the same seed while also still acting random wrt time.
	x = sim.run(n_runs=1, method='ode')
	dtf = x.dataframe
	cp_params = np.zeros_like(new_params[cell_ind])
	cp_params[:] = new_params[cell_ind, :]
	cp_params[12] = dtf["obs_IGFR"].iloc[-1]
	cp_params[14] = dtf["obs_FoxO_c"].iloc[-1]
	cp_params[15] = dtf["obs_FoxO_n"].iloc[-1]
	fig,ax = plt.subplots(1)
	for l in range(len(L)):
		cp_params[16] = L[l]
		sim = BngSimulator(model, tspan=t_main, param_values=cp_params)
		x = sim.run(n_runs=1, method='ode')
		dtf = x.dataframe
		ax.plot(t_main/60,dtf["obs_FoxO_n"],color = dose_colors[l])
		sim = BngSimulator(model, tspan=t_main, param_values=cp_params)
		x = sim.run(n_runs=1, method='ssa')
		dtf = x.dataframe
		ax.plot(t_main / 60, dtf["obs_FoxO_n"],color = dose_colors[l],linestyle = 'dashed')
	ax.spines["top"].set_visible(False)
	ax.spines["right"].set_visible(False)
	ax.set_ylabel("Nuclear Foxo")
	ax.set_xlabel("Time (Minutes)")
	plt.tight_layout()
	plt.show()

