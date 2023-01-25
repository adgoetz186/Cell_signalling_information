import numpy as np
import time
import scipy.stats as st
from Mutual_Information_Final_Version.functions.toy_model_functions import Toy_Model_Functions as tmf
import Mutual_Information_Final_Version.functions.file_processing_functions.user_input_functions as uif
import pandas as pd

# Simulates population of cells with networks defined by toy receptor model for different levels of parameter variation
# Records:
# Single cell mutual information for population of cells for various distributions of parameters
# Mutual information of cell state agnostic response for various distributions of parameters

# _____ Load arguments BEGIN _____
# Specifies path to file which specifies user specified arguments
user_input_file_path = "/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Figure_Generating_Programs/Figure_2/General_Data_Generation/Generate_Framework_Comparison_Toy_Model/args.txt"
user_arguments = uif.read_file_to_arg_dict(user_input_file_path)
# _____ Load arguments END _____


# _____ File path declarations BEGIN _____
# location of the output data
output_file_path = user_arguments["output_file_path"]
# output data filename prefix
output_file_name_prefix = user_arguments["output_file_name_prefix"]
ceemi_deg_array_filename = f"{output_file_path}/{output_file_name_prefix}_cmi_deg.csv"
ceemi_r0_array_filename = f"{output_file_path}/{output_file_name_prefix}_cmi_r0.csv"
mi_deg_array_filename = f"{output_file_path}/{output_file_name_prefix}_mi_deg.csv"
mi_r0_array_filename = f"{output_file_path}/{output_file_name_prefix}_mi_r0.csv"
# _____ File path declarations END _____

# _____ Main Code BEGIN _____
# cv of input distribution
input_cv = user_arguments["input_cv"]

# mean of input distribution
input_mean = user_arguments["input_mean"]

# Number of bins to use for input discretization
input_partitions = user_arguments["input_partitions"]

# Values of parameter cv's to plot, linspace arguement, base 10
param_cv_linspace_arg = user_arguments["param_cv_linspace_arg"]
parameter_cv = 10**np.linspace(param_cv_linspace_arg[0],param_cv_linspace_arg[1],param_cv_linspace_arg[2])

# average degredation rate, k_deg
deg_mean = user_arguments["deg_mean"]

# average initial receptor count value, r0
r0_mean = user_arguments["r0_mean"]

# rate of binding, k_bind
k_bind = user_arguments["k_bind"]

# rate of unbinding, k_unbind
k_unbind = user_arguments["k_unbind"]

# for case where degradation varies
r0_val = user_arguments["r0_val"]

# for case where receptor average varies
deg_val = user_arguments["deg_val"]

# Performs this many monte carlo draws for each integral
mc_draws = user_arguments["mc_draws"]

# MI array for different values of r0 variance
mir0 = np.zeros((1,np.size(parameter_cv)))

# Cee-MI array for different values of r0 variance
cmi_r0_array = np.zeros((mc_draws,np.size(parameter_cv)))

# MI array for different values of kdeg variance
mideg = np.zeros((1,np.size(parameter_cv)))

# Cee-MI array for different values of kdeg variance
cmi_deg_array = np.zeros((mc_draws,np.size(parameter_cv)))

# obtains the input variance
input_variance = (input_cv*input_mean)**2

# Generates input distribution by percentile binning of a gamma
u_scale1 = input_variance/ input_mean
uvar_shape1 = input_mean / u_scale1
edge = np.linspace(0, 1, input_partitions+1)
center = (edge[1:]+edge[:-1])/2
ulist = st.gamma.ppf(center, uvar_shape1, scale=u_scale1)

# defines the kdeg variance from the mean and CV
deg_variance = (deg_mean * parameter_cv) ** 2

# defines the r0 variance from the mean and CV
r0_variance = (r0_mean * parameter_cv) ** 2

# used to track estimated time remaining
time_start = time.time()

for pvarind in range(np.size(parameter_cv)):
    # Calculates the Cee-MI for the given kdeg variance
    cmi_deg_array[:,pvarind] = -1*tmf.cmi_deg_var(ulist,deg_variance[pvarind],deg_mean,mc_draws,r0_val,k_bind,k_unbind)

    # Calculates the Cee-MI for the given r0 variance
    cmi_r0_array[:,pvarind] = -1*tmf.cmi_r0_var(ulist,r0_variance[pvarind],r0_mean,mc_draws,deg_val,k_bind,k_unbind)

    # Calculates the MI of the CSAR for the given deg variance
    mideg[0,pvarind] = -1*tmf.mi_deg_var(ulist,deg_variance[pvarind],deg_mean,mc_draws,r0_val,k_bind,k_unbind)

    # Calculates the MI of the CSAR for the given r0 variance
    mir0[0,pvarind] = -1*tmf.mi_r0_var(ulist,r0_variance[pvarind],r0_mean,mc_draws,deg_val,k_bind,k_unbind)
    
    # Tracks time remaining
    time_elapsed = time.time() - time_start
    print(f"Approximately {round(time_elapsed/(pvarind+1)*((np.size(input_variance)-1)*np.size(parameter_cv)+np.size(parameter_cv)-pvarind-1)/60)} minutes remain")
    
# Saves resulting arrays
pd.DataFrame(cmi_deg_array,columns=parameter_cv).to_csv(ceemi_deg_array_filename)
pd.DataFrame(cmi_r0_array,columns=parameter_cv).to_csv(ceemi_r0_array_filename)
pd.DataFrame(mideg,columns=parameter_cv).to_csv(mi_deg_array_filename)
pd.DataFrame(mir0,columns=parameter_cv).to_csv(mi_r0_array_filename)
# _____ Main Code END _____
