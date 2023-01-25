import numpy as np
import scipy.stats as st
from Mutual_Information_Final_Version.functions.toy_model_functions import Toy_Model_Functions as tmf
import Mutual_Information_Final_Version.functions.file_processing_functions.user_input_functions as uif
import pandas as pd

# Simulates single cell toy receptor model for different cell parameter values
# Records:
# Single cell mutual information for various cell parameter values

# _____ Load arguments BEGIN _____
# Specifies path to file which specifies user specified arguments
user_input_file_path = "/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Figure_Generating_Programs/Figure_2/General_Data_Generation/Generate_Single_Cell_Performances/args.txt"
user_arguments = uif.read_file_to_arg_dict(user_input_file_path)
# _____ Load arguments END _____

# _____ File path declarations BEGIN _____
output_file_path = user_arguments["output_file_path"]
output_file_name_prefix = user_arguments["output_file_name_prefix"]
output_dir_deg = f"{output_file_path}/{output_file_name_prefix}_cmi_single_cell_deg.csv"
output_dir_r0 = f"{output_file_path}/{output_file_name_prefix}_cmi_single_cell_r0.csv"
# _____ File path declarations END _____


# _____ Figure generation BEGIN _____
input_cv = user_arguments["input_cv"]

input_mean = user_arguments["input_mean"]
input_partitions = user_arguments["input_partitions"]

param_cv_linspace_arg = user_arguments["param_cv_linspace_arg"]
param_cv = 10 ** np.linspace(param_cv_linspace_arg[0], param_cv_linspace_arg[1], param_cv_linspace_arg[2])
# parameter_cv = np.linspace(-1,1,15)
deg_mean = user_arguments["deg_mean"]
r0_mean = user_arguments["r0_mean"]
k_bind = user_arguments["k_bind"]
k_unbind = user_arguments["k_unbind"]

# for case where degradation varies
r0_val = user_arguments["r0_val"]

# for case where receptor average varies
deg_val = user_arguments["deg_val"]

# Performs this many monte carlo draws for each integral
mc_draws = user_arguments["mc_draws"]

# The number of parameter points to evaluate the single cell performance
parameter_partitions = user_arguments["parameter_partitions"]

# range of k_deg
kdr = user_arguments["kdr"]

# range of r0
r0r = user_arguments["r0r"]

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
