from Mutual_Information_Main.functions.Models.Toy import model
import numpy as np
import os
from pathlib import Path, PureWindowsPath
from pysb.simulator import BngSimulator
import multiprocessing as mp
import functools
import time
import random
import scipy.stats as st


# This program generates outcomes of the EGFR pathway using the Gillespie algorithm
# These results are used to test the moment closure accuracy

def single_run(r, params, t_main, L):
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
	total_cell = np.shape(params)[0]
	cp_params = np.zeros_like(params)
	cp_params[:] = params
	# The initial receptor count is the result of a steady state creation annihilation event and is poisson as a result
	cp_params[4] = st.poisson.rvs(cp_params[0]/cp_params[1])
	for l in range(len(L)):
		cp_params[5] = L[l]
		sim = BngSimulator(model, tspan=t_main, param_values=cp_params)
		# The seed is incredibly important for parallel processing since the time is used to create random seeds
		# often all parallel workers will end up with the same time which makes them all give identical results
		# this prevents workers from having the same seed while also still acting random wrt time.
		x = sim.run(n_runs=1, method='ssa',seed = r+random.randint(0,int(2**63-1-r-total_cell)))
		dtf = x.dataframe
		cell_run_array.append(dtf["obs_B"].iloc[-1])
	return cell_run_array


if __name__ == '__main__':
	# _____ Setting the CWD to be Mutual_Information_Main BEGIN _____
	# Cell_signaling_information path here
	path_to_CSI = ""
	if path_to_CSI == "":
		try:
			# Obtains the location of the Cell_signaling_information folder if it is in the cwd parents
			path_to_CSI = Path.cwd().parents[
				[Path.cwd().parents[i].parts[-1] for i in range(len(Path.cwd().parts) - 1)].index(
					"Cell_signalling_information")]
		except ValueError:
			print("Cell_signalling_information not found in cwd parents, trying sys.path")
			try:
				# Obtains the location of the Cell_signaling_information folder if it is in sys.path
				path_to_CSI = Path(sys.path[[Path(i).parts[-1] for i in sys.path].index("Cell_signalling_information")])
			except ValueError:
				print("Cell_signalling_information not found in sys.path "
				      "consult 'Errors with setting working directory' in README")
	else:
		path_to_CSI = Path(path_to_CSI)
	path_to_MIM = path_to_CSI / "Mutual_Information_Main"
	os.chdir(path_to_MIM)
	# _____ Setting the CWD to be Mutual_Information_Main END _____
	
	# _____ File path declarations BEGIN _____
	# Specifies path to output directory
	output_folder = Path("Data/Toy_Model/Gillespie_Results/")
	# _____ File path declarations END _____
	
	
	# _____ Main code BEGIN _____
	
	# average degredation rate, k_deg
	deg_mean = 5
	# average initial receptor count value, r0
	r0_mean = 500
	# rate of binding, k_bind
	k_bind = 1
	# rate of unbinding, k_unbind
	k_unbind = 10
	
	input_mean = 10
	input_cv = 1
	# obtains the input variance
	input_variance = (input_cv * input_mean) ** 2
	input_partitions = 25
	# Generates input distribution by percentile binning of a gamma
	u_scale1 = input_variance / input_mean
	uvar_shape1 = input_mean / u_scale1
	edge = np.linspace(0, 1, input_partitions + 1)
	center = (edge[1:] + edge[:-1]) / 2
	ulist = st.gamma.ppf(center, uvar_shape1, scale=u_scale1)
	
	# time of main simulation
	t_f = 3600
	# recorded time steps
	t_main = np.arange(0, t_f + t_f, t_f)
	# number of simulations to run per cell (2500)
	n_sim = 2500

	# number of cells to use
	cell_count = 1
	#pars = [ksyn, k1, kn1, kap, kdp, ki, kis]

	k_prod = deg_mean*r0_mean
	params = np.array([k_prod,deg_mean,k_bind,k_unbind,0.0,0.0])
	time_start = time.time()
	# Runs the ssa algorithm
	for cell_ind in range(cell_count):
		# uses all but 1 cpu, this is done to allow 1 cpu for other programs
		with mp.Pool(mp.cpu_count() - 1) as p:
			pfun = functools.partial(single_run, params=params,
			                         t_main=t_main, L=ulist)
			cell_array = np.array(p.map(pfun, range(n_sim)))
		np.save(output_folder / f"cell_{cell_ind}", cell_array)
		print((time.time() - time_start) / (60 * (cell_ind + 1)))
	# generates header dictionaries
	cell_header = []
	for l in range(len(ulist)):
		cell_header.append({'dose': ulist[l]})
	with open(output_folder / "cell_header.txt", 'w') as co:
		for i in cell_header:
			co.write(str(i) + '\n')
# _____ Main code END _____
