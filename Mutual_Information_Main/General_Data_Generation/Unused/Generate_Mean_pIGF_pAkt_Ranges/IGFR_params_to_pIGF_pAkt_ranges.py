from Mutual_Information_Final_Version.functions.Models.IGFR_Akt_FoxO import model
import numpy as np
from pysb.simulator import BngSimulator
import multiprocessing as mp
import functools
import Mutual_Information_Final_Version.functions.file_processing_functions.user_input_functions as uif
import time
import random
import copy as cp

# This program generates outcomes of the IGFR/Akt/FoxO pathway using the Gillespie algorithm
# These results are used to test the moment closure accuracy
	


# _____ Load arguments BEGIN _____
# Specifies path to file which specifies user specified arguments, args file should be in same folder as this program
user_input_file_path = "args.txt"
user_arguments = uif.read_file_to_arg_dict(user_input_file_path)
# _____ Load arguments END _____

# _____ File path declarations BEGIN _____
# Specifies path to drawn model cell parameters
params_file_path = user_arguments["params_file_path"]
# Specifies path to output directory
output_dir_path = user_arguments["output_dir_path"]
# _____ File path declarations END _____

# _____ Loading files BEGIN _____
# Loads parameters of drawn cells
params = 10**np.loadtxt(params_file_path, delimiter=",")
# _____ Loading files END _____

# _____ Main code BEGIN _____
# time allowed before sampling cells with no dose
times = np.array(user_arguments["times"])*60

# list of doses
L = user_arguments["L"]
# number of cells to use
cell_count = user_arguments["cell_count"]



# The set of parameters used elsewhere doesnt allow for non equilibrium initial conditions
# this transforms the parameters into a set which can allow for a relevant subset of
# non equilibrium initial conditions
kprod = params[:,0]*params[:,1]
akt_0 = params[:,12]
n_foxo_0 = (params[:,13]*params[:,10]/(params[:,10]+params[:,11]))
c_foxo_0 = (params[:,13]*params[:,11]/(params[:,10]+params[:,11]))
igfr_0 = params[:,0]
new_params = np.zeros((np.shape(params)[0],np.shape(params)[1]+3))
new_params[:,0] = params[:,0]*params[:,1]
new_params[:,1:12] = params[:,1:12]
new_params[:,12] = igfr_0
new_params[:,13] = akt_0
new_params[:,14] = c_foxo_0
new_params[:,15] = n_foxo_0
pIGF_array = np.zeros((cell_count,np.size(times)*len(L)))
pAkt_array = np.zeros((cell_count, np.size(times)*len(L)))
nfoxo_array = np.zeros((cell_count, np.size(times)*len(L)))
# Runs the ssa algorithm
s = time.time()
for cell_ind in range(cell_count):
	print(cell_ind)
	pIGF = []
	pAkt = []
	nfoxo = []
	for l in range(len(L)):
		cp_params = np.zeros_like(new_params[cell_ind])
		cp_params[:] = new_params[cell_ind, :]
		cp_params[16] = L[l]
		sim = BngSimulator(model, tspan=times, param_values=cp_params)
		x = sim.run(n_runs=1, method='ode')
		dtf = x.dataframe
		for t_ind in range(np.size(times)):
			pIGF.append(dtf["obs_pIGFR"].iloc[t_ind])
			pAkt.append(dtf["obs_pAkt"].iloc[t_ind])
			nfoxo.append(dtf["obs_FoxO_n"].iloc[t_ind])
	pIGF_array[cell_ind] = np.array(pIGF)
	pAkt_array[cell_ind] = np.array(pAkt)
	nfoxo_array[cell_ind] = np.array(nfoxo)
print(time.time() - s)
np.savetxt(f"{output_dir_path}pIGF.csv", pIGF_array, delimiter=',')
np.savetxt(f"{output_dir_path}pAkt.csv", pAkt_array, delimiter=',')
np.savetxt(f"{output_dir_path}nFoxo.csv", nfoxo_array, delimiter=',')
# generates header dictionaries
cell_header = []
for l in range(len(L)):
	for t_ind in range(np.size(times)):
		cell_header.append({"time": times[t_ind], 'dose': L[l]})
with open(f"{output_dir_path}/cell_header.txt",'w') as co:
	for i in cell_header:
		co.write(str(i)+'\n')
# _____ Main code END _____
