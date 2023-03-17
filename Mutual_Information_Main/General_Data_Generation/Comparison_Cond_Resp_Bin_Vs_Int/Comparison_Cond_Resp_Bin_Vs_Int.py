import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot
import sys
import os
from pathlib import Path, PureWindowsPath
import Mutual_Information_Main.functions.mutual_information_functions.MI_Calculation_Functions as mi
import time
# Calculates single cell mutual information values using the binning method used elsewhere in this project and
# a more direct numerical integration using scipy.integrate.quad. This is done to ensure binning method accuracy.

# _____ Setting the CWD to be Mutual_Information_Main BEGIN _____
# Cell_signaling_information path here
path_to_CSI = ""
if path_to_CSI == "":
	try:
		# Obtains the location of the Cell_signaling_information folder if it is in the cwd parents
		path_to_CSI = Path.cwd().parents[[Path.cwd().parents[i].parts[-1] for i in range(len(Path.cwd().parts) - 1)].index(
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
path_to_MIM = path_to_CSI/"Mutual_Information_Main"
os.chdir(path_to_MIM)
# _____ Setting the CWD to be Mutual_Information_Main END _____


# _____ File path declarations BEGIN _____
# location of the output egfr data
output_egfr_files = Path("Data/EGFR/moment_comparison_integration_methods/")

# location of the output igfr data
output_igfr_files = Path("Data/IGFR/moment_comparison_integration_methods/")

# location of egfr moment doses
egfr_dose_moment_path = Path("Data/EGFR/Moments/noise_factor_1/moments.csv")

# location of igfr moment doses
igfr_dose_moment_path = Path("Data/IGFR/Moments/Model_Moments/4_dose_response/moments_162_4_dose_60_90_min.csv")
# _____ File path declarations END _____

# _____ Loading files BEGIN _____
# Loads egfr and igfr dose response moments
egfr_dose_moment = np.loadtxt(egfr_dose_moment_path,delimiter=',')
igfr_dose_moment = np.loadtxt(igfr_dose_moment_path,delimiter=',')
# _____ Loading files END _____

# _____ Data Processing BEGIN _____
relative_input = np.ones(np.shape(egfr_dose_moment)[1]//2)
relative_input/=np.size(relative_input)
val = np.zeros((1,np.size(relative_input)*2))

number_to_run = 50
rand_int = np.arange(np.shape(egfr_dose_moment)[0])
np.random.shuffle(rand_int)
rand_int = rand_int[:number_to_run]
moment_array_to_use = egfr_dose_moment[rand_int]

# Creates and saves population level information calculation for both the binned and integration approaches
population_average_moments = np.average(moment_array_to_use,axis=0)
pop_int = np.ones(1)*mi.mutual_information_direct_integration_gamma(relative_input, np.reshape(population_average_moments,(1,-1)))[0]

pop_bin = np.ones(1)*-1*mi.conditional_mutual_information_from_list_of_moments(relative_input, [np.reshape(population_average_moments,(1,-1))], assumed_distribution="gamma", discretization_parameter=0.05, percentile_cutoff=0.0005)
np.save(output_egfr_files/"pop_integration",pop_int)
np.save(output_egfr_files/"pop_binning",pop_bin)

egfr_times_array = np.zeros(2)

integral_array = np.zeros(number_to_run)
binned_array = np.zeros(number_to_run)

start = time.time()
for i in range(np.shape(moment_array_to_use)[0]):
	print(1,i)
	val[0] = moment_array_to_use[i]
	binned_array[i] = -1*mi.conditional_mutual_information_from_list_of_moments(relative_input, [val], assumed_distribution="gamma", discretization_parameter=0.05, percentile_cutoff=0.0005)
egfr_times_array[0] = time.time()-start
start = time.time()
for i in range(np.shape(moment_array_to_use)[0]):
	print(2,i)
	val[0] = moment_array_to_use[i]
	integral_array[i] = mi.mutual_information_direct_integration_gamma(relative_input, val)[0]
	# scipy.integrate.quad can struggle with small pulses and because of this can miss a conditional response,
	# as a result if the values of both methods are sufficiently different a more stringent error term is used for
	# the quad function to attempt to correct this error the process is slwoer so is only invoked when needed.
	# This problem is only significant for single cell egfr as the responses are so narrow wrt the integration range
	if abs(integral_array[i]-binned_array[i])/integral_array[i]*100 >1:
		integral_array[i] = mi.mutual_information_direct_integration_gamma(relative_input, val,error=1e-30)[0]
		if abs(integral_array[i] - binned_array[i]) / integral_array[i] * 100 > 1:
			print(":(")
			input()
egfr_times_array[1] = time.time()-start
np.save(output_egfr_files/"sc_integration",integral_array)
np.save(output_egfr_files/"sc_binning",binned_array)
np.save(output_egfr_files/"time_of_runs",egfr_times_array)







igfr_avg_moment = (igfr_dose_moment[:,0:-1:2]+igfr_dose_moment[:,1::2])/2
relative_input = np.ones(4)/4
rand_int = np.arange(np.shape(igfr_avg_moment)[0])
np.random.shuffle(rand_int)
rand_int = rand_int[:number_to_run]
moment_array_to_use = igfr_avg_moment[rand_int]
# Creates and saves population level information calculation for both the binned and integration approaches
population_average_moments = np.average(moment_array_to_use,axis=0)
pop_int = np.ones(1)*mi.mutual_information_direct_integration_gamma(relative_input, np.reshape(population_average_moments,(1,-1)))[0]
pop_bin = np.ones(1)*-1*mi.conditional_mutual_information_from_list_of_moments(relative_input, [np.reshape(population_average_moments,(1,-1))], assumed_distribution="gamma", discretization_parameter=0.05, percentile_cutoff=0.0005)
np.save(output_igfr_files/"pop_integration",pop_int)
np.save(output_igfr_files/"pop_binning",pop_bin)

igfr_times_array = np.zeros(2)

val = np.zeros((1,np.size(relative_input)*2))
start = time.time()

binned_array = np.zeros(number_to_run)
start = time.time()
for i in range(number_to_run):
	print(3,i)
	val[0] = igfr_avg_moment[i]
	binned_array[i] = -1*mi.conditional_mutual_information_from_list_of_moments(relative_input, [val], assumed_distribution="gamma", discretization_parameter=0.05, percentile_cutoff=0.0005)
igfr_times_array[0] = time.time()-start

integral_array = np.zeros(number_to_run)
start = time.time()
for i in range(number_to_run):
	print(4,i)
	val[0] = igfr_avg_moment[i]
	integral_array[i] = mi.mutual_information_direct_integration_gamma(relative_input, val)[0]
igfr_times_array[1] = time.time()-start

np.save(output_igfr_files/"sc_integration",integral_array)
np.save(output_igfr_files/"sc_binning",binned_array)
np.save(output_igfr_files/"time_of_runs",igfr_times_array)
# _____ Data Processing END _____
