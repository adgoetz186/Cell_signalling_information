import numpy as np
import time
import os
import scipy.stats as st
from pathlib import Path, PureWindowsPath
from Mutual_Information_Main.functions.toy_model_functions import Toy_Model_Functions as tmf
import pandas as pd

# Simulates population of cells with networks defined by toy receptor model for different levels of parameter variation
# Records:
# Single cell mutual information for population of cells for various distributions of parameters
# Mutual information of cell state agnostic response for various distributions of parameters

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
# location of the output data
output_file_path = Path("Data/Toy_Model/Information_Values/framework_comparison")
# output data filename prefix
output_file_name_prefix = "toy_model_system"
ceemi_deg_array_filename = output_file_path/(output_file_name_prefix+"_cmi_deg.csv")
ceemi_r0_array_filename = output_file_path/(output_file_name_prefix+"_cmi_r0.csv")
mi_deg_array_filename = output_file_path/(output_file_name_prefix+"_mi_deg.csv")
mi_r0_array_filename = output_file_path/(output_file_name_prefix+"_mi_r0.csv")
# _____ File path declarations END _____

# _____ Main Code BEGIN _____
# cv of input distribution
input_cv = 1

# mean of input distribution
input_mean = 10

# Number of bins to use for input discretization
input_partitions = 25

# Values of parameter cv's to plot, linspace arguement, base 10
parameter_cv = 10**np.linspace(-1.5,-0.5,50)

# average degredation rate, k_deg
deg_mean = 5

# average initial receptor count value, r0
r0_mean = 500

# rate of binding, k_bind
k_bind = 1

# rate of unbinding, k_unbind
k_unbind = 10

# for case where degradation varies
r0_val = 50

# for case where receptor average varies
deg_val = 5

# Performs this many monte carlo draws for each integral (50000)
mc_draws = 50000

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
