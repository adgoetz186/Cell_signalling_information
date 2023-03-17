import numpy as np
import scipy.stats as st
from pathlib import Path, PureWindowsPath
import os
from Mutual_Information_Main.functions.toy_model_functions import Toy_Model_Functions as tmf
import pandas as pd

# Simulates single cell toy receptor model for different cell parameter values
# Records:
# Single cell mutual information for various cell parameter values

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
output_file_path = Path("Data/Toy_Model/Information_Values/single_cell_performance")
output_file_name_prefix = "single_cell_toy_model_system"
output_dir_deg = output_file_path/(output_file_name_prefix + "_cmi_single_cell_deg.csv")
output_dir_r0 = output_file_path/(output_file_name_prefix + "_cmi_single_cell_r0.csv")
# _____ File path declarations END _____


# _____ Figure generation BEGIN _____
input_cv = 1

input_mean = 10
input_partitions = 25

deg_mean = 5
r0_mean = 500
k_bind = 1
k_unbind = 10

# for case where degradation varies
r0_val = 50

# for case where receptor average varies
deg_val = 5

# Performs this many monte carlo draws for each integral (50000)
mc_draws = 50000

# The number of parameter points to evaluate the single cell performance
parameter_partitions = 100

# range of k_deg
kdr = [0,2]

# range of r0
r0r = [0,3.5]

# creates parameter arrays
log_deg_mean = np.linspace(kdr[0], kdr[1], parameter_partitions)
log_r0_mean = np.linspace(r0r[0], r0r[1], parameter_partitions)

# Creates arrays to store results
single_cell_performances_deg = np.zeros((1, parameter_partitions))
single_cell_performances_r0 = np.zeros((1, parameter_partitions))

# Creates array of input variances
input_variance = (input_cv * input_mean) ** 2

# obtains input distribution
u_scale1 = input_variance / input_mean
uvar_shape1 = input_mean / u_scale1
edge = np.linspace(0, 1, input_partitions + 1)
center = (edge[1:] + edge[:-1]) / 2
ulist = st.gamma.ppf(center, uvar_shape1, scale=u_scale1)

# obtains single cell performances
for parm_ind in range(parameter_partitions):
	single_cell_performances_deg[0,parm_ind] = tmf.single_cell_MI(ulist, 10**log_deg_mean[parm_ind], r0_val, k_bind, k_unbind)
	single_cell_performances_r0[0, parm_ind] = tmf.single_cell_MI(ulist, deg_val, 10 ** log_r0_mean[parm_ind], k_bind,k_unbind)

# saves results
pd.DataFrame(single_cell_performances_deg, index=[input_cv], columns=log_deg_mean).to_csv(output_dir_deg)
pd.DataFrame(single_cell_performances_r0, index=[input_cv], columns=log_r0_mean).to_csv(output_dir_r0)
# _____ Figure generation END _____
