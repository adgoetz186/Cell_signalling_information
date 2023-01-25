from Mutual_Information_Final_Version.functions.Models.EGFR import model
import numpy as np
from pysb.simulator import BngSimulator
import multiprocessing as mp
import functools
import Mutual_Information_Final_Version.functions.file_processing_functions.user_input_functions as uif
import time
import random
import scipy.stats as st


# This program generates outcomes of the EGFR pathway using the Gillespie algorithm
# These results are used to test the moment closure accuracy

def single_run(r, new_params, cell_ind, t_main, L):
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
	print(r)
	cell_run_array = []
	total_cell = np.shape(new_params)[0]
	cp_params = np.zeros_like(new_params[cell_ind])
	cp_params[:] = new_params[cell_ind, :]
	# The initial receptor count is the result of a steady state creation annihilation event and is poisson as a result
	cp_params[7] = st.poisson.rvs(cp_params[0]/cp_params[5])
	for l in range(len(L)):
		cp_params[8] = L[l]
		sim = BngSimulator(model, tspan=t_main, param_values=cp_params)
		# The seed is incredibly important for parallel processing since the time is used to create random seeds
		# often all parallel workers will end up with the same time which makes them all give identical results
		# this prevents workers from having the same seed while also still acting random wrt time.
		x = sim.run(n_runs=1, method='ssa',seed = r+random.randint(0,int(2**63-1-r-total_cell)))
		dtf = x.dataframe
		cell_run_array.append(dtf["obs_sEGFR"].iloc[-1])
	return cell_run_array


if __name__ == '__main__':
	# _____ Load arguments BEGIN _____
	# Specifies path to file which specifies user specified arguments
	user_input_file_path = "/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Figure_Generating_Programs/S_Gillespie_Supplimentary/General_Data_Generation/EGFR/ssa_EGFR_args.txt"
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
	# time of main simulation
	t_f = user_arguments["t_f"]
	# recorded time steps
	t_main = np.arange(0, t_f + t_f, t_f)
	# number of simulations to run per cell
	n_sim = user_arguments["n_sim"]
	# list of doses
	L = user_arguments["L"]
	# number of cells to use
	cell_count = user_arguments["cell_count"]
	#pars = [ksyn, k1, kn1, kap, kdp, ki, kis]
	
	# scales synthesis rate
	params[:, 0] /= 0.00122
	
	EGFR_0 = (params[:, 0] / params[:, 5]).astype(int)
	new_params = np.zeros((np.shape(params)[0], np.shape(params)[1] + 2))
	new_params[:,:7] = params
	new_params[:, 7] = EGFR_0
	time_start = time.time()
	# Runs the ssa algorithm
	for cell_ind in range(cell_count):
		# uses all but 1 cpu, this is done to allow 1 cpu for other programs
		with mp.Pool(mp.cpu_count() - 1) as p:
			pfun = functools.partial(single_run, new_params=new_params, cell_ind=cell_ind,
			                         t_main=t_main, L=L)
			cell_array = np.array(p.map(pfun, range(n_sim)))
		np.savetxt(f"{output_folder}/cell_{cell_ind}.csv", cell_array, delimiter=',')
		print((time.time() - time_start) / (60 * (cell_ind + 1)))
	# generates header dictionaries
	cell_header = []
	for l in range(len(L)):
		cell_header.append({'dose': L[l]})
	with open(f"{output_folder}/cell_header.txt", 'w') as co:
		for i in cell_header:
			co.write(str(i) + '\n')
# _____ Main code END _____
