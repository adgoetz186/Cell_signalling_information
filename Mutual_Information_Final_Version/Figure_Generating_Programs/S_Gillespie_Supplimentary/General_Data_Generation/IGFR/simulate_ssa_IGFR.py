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

def single_run(r,t_f1_array,new_params,cell_ind,t_main,L):
	"""performs single run of the system involving 2 coupled Gillespie algorithms
	
	Input:
	
		-r (int): the run number
		
		-t_f1_array (ndarray): the array containing the duration of the initial simulation
		
		-new_params (ndarray): the array containing all cell parameters, with each cell having a unique row
		
		-cell_ind (int): specifies the cell currently being simulates, indexes the row of the param array
		
		-t_main (ndarray): specifies all times to generate and store results for the main simulation
		
		-L (list): list of all doses UNITS ARE IMPORTANT and are assumed to be nM
		
	Output:
	
		-cell_run_array (list): a single trajectory of nuclear foxo values"""
	cell_run_array = []
	total_cell = np.shape(new_params)[0]
	sim = BngSimulator(model, tspan=t_f1_array, param_values=new_params[cell_ind])
	# The seed is incredibly important for parallel processing since the time is used to create random seeds
	# often all parallel workers will end up with the same time which makes them all give identical results
	# this prevents workers from having the same seed while also still acting random wrt time.
	x = sim.run(n_runs=1, method='ssa',seed = r+random.randint(0,int(2**63-1-r)))
	dtf = x.dataframe
	cp_params = np.zeros_like(new_params[cell_ind])
	cp_params[:] = new_params[cell_ind,:]
	cp_params[12] = dtf["obs_IGFR"].iloc[-1]
	cp_params[14] = dtf["obs_FoxO_c"].iloc[-1]
	cp_params[15] = dtf["obs_FoxO_n"].iloc[-1]
	for l in range(len(L)):
		cp_params[16] = L[l]
		sim = BngSimulator(model, tspan=t_main, param_values=cp_params)
		x = sim.run(n_runs=1, method='ssa',seed = r+total_cell+random.randint(0,int(2**63-1-r-total_cell)))
		dtf = x.dataframe
		for t_ind in range(np.size(t_main)):
			cell_run_array.append(dtf["obs_FoxO_n"].iloc[t_ind])
	return cell_run_array
	

if __name__ == '__main__':
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
	params = 10**np.loadtxt(params_file_path, delimiter=",")
	# _____ Loading files END _____
	
	# _____ Main code BEGIN _____
	# time allowed before sampling cells with no dose
	t_f1 = user_arguments["t_0"]
	t_f1_array = np.linspace(0,t_f1,2)
	
	# time of main simulation
	t_f2 = user_arguments["t_f2"]
	# recorded time steps
	t_delta = user_arguments["t_delta"]
	t_main = np.arange(0,t_f2+t_delta,t_delta)
	# number of simulations to run per cell
	n_sim = user_arguments["n_sim"]
	# list of doses
	L = user_arguments["L"]
	# number of cells to use
	cell_count = user_arguments["cell_count"]

	
	
	# The set of parameters used elsewhere doesnt allow for non equilibrium initial conditions
	# this transforms the parameters into a set which can allow for a relevant subset of
	# non equilibrium initial conditions
	kprod = params[:,0]*params[:,1]
	akt_0 = params[:,12].astype(int)
	n_foxo_0 = (params[:,13]*params[:,10]/(params[:,10]+params[:,11])).astype(int)
	c_foxo_0 = (params[:,13]*params[:,11]/(params[:,10]+params[:,11])).astype(int)+1
	igfr_0 = params[:,0].astype(int)
	new_params = np.zeros((np.shape(params)[0],np.shape(params)[1]+3))
	new_params[:,0] = params[:,0]*params[:,1]
	new_params[:,1:12] = params[:,1:12]
	new_params[:,12] = igfr_0
	new_params[:,13] = akt_0
	new_params[:,14] = c_foxo_0
	new_params[:,15] = n_foxo_0
	
	# Runs the ssa algorithm
	for cell_ind in range(cell_count):
		# uses all but 1 cpu, this is done to allow 1 cpu for other programs
		with mp.Pool(mp.cpu_count()-1) as p:
			pfun = functools.partial(single_run,t_f1_array = t_f1_array,new_params = new_params,cell_ind = cell_ind,t_main = t_main,L = L)
			cell_array = p.map(pfun, range(n_sim))
		np.savetxt(f"{output_folder}/cell_{cell_ind}.csv",cell_array,delimiter=',')
	# generates header dictionaries
	cell_header = []
	for l in range(len(L)):
		for t_ind in range(np.size(t_main)):
			cell_header.append({"time": t_main[t_ind], 'dose': L[l]})
	with open(f"{output_folder}/cell_header.txt",'w') as co:
		for i in cell_header:
			co.write(str(i)+'\n')
	# _____ Main code END _____
