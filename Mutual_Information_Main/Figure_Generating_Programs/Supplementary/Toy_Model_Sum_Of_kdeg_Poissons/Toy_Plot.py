import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
from Mutual_Information_Main.functions.toy_model_functions import Toy_Model_Functions as tmf
from matplotlib.patches import Rectangle
from pathlib import Path
import os

# Generates plot of EGFR signaling performance with difference levels of noise and MI of cell state agnostic response

# sets font for figure text
plt.rcParams.update({'font.size': 15})

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
# location of output image
output_path = Path("Figures/Supplementary/sum_of_toy_responses_kdeg/")

# name of output image
output_name = "sum_of_toy_responses_kdeg"
# _____ File path declarations END _____

# _____ Main Code BEGIN _____
# dimensions of figure
fig_dim = [4.2,8]
# cv of input distribution
input_cv = 1

# mean of input distribution
input_mean = 10

# Number of bins to use for input discretization
input_partitions = 25

# Values of parameter cv
parameter_cv = 10**-.5

# average degredation rate, k_deg
deg_mean = 5

# average initial receptor count value, r0
r0_val = 50

# rate of binding, k_bind
k_bind = 1

# rate of unbinding, k_unbind
k_unbind = 10

# Performs this many monte carlo draws for each integral (50000)
mc_draws = 50000

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

# generate population response under assumption of single cell poisson and then average
sc_crm_list = tmf.pxgut_deg(ulist,deg_variance,deg_mean,mc_draws,r0_val,k_bind,k_unbind,5)
xlist = np.arange(r0_val*5)
direct_average_poisson = np.zeros_like(sc_crm_list[0])
count = 0
for i in sc_crm_list:
	direct_average_poisson += i
direct_average_poisson /= len(sc_crm_list)

mean_sc_response = np.zeros((np.size(ulist),len(sc_crm_list)))
for i in range(len(sc_crm_list)):
	# Calculation of sc means and vars from crms
	# The single cell responses are defined as poisson so mean = variance
	mean_sc_response[:,i] = np.sum(xlist*sc_crm_list[i],axis=1)


average_mean = np.average(mean_sc_response,axis=1)
average_var = np.average(mean_sc_response + mean_sc_response**2,axis=1) - average_mean**2

fig,axs = plt.subplots(9,1,sharex = True, figsize=(fig_dim[0], fig_dim[1]))
for subplot_number in range(9):
	i = subplot_number*3
	p = average_mean[i]/average_var[i]
	n = average_mean[i]**2/(average_var[i]-average_mean[i])
	#shift = np.sum(xlist*direct_average_poisson[i])
	pxgutdist = st.nbinom.pmf(xlist,n,p)
	axs[subplot_number].plot(xlist,direct_average_poisson[i],label= "Avg. Poisson", color = 'k')
	axs[subplot_number].plot(xlist,pxgutdist, label = "nBinom", color = "red",linestyle = "--")
	axs[subplot_number].spines["top"].set_visible(False)
	axs[subplot_number].spines["right"].set_visible(False)
	axs[subplot_number].spines["left"].set_visible(False)
	axs[subplot_number].axes.yaxis.set_visible(False)
	label = f"{round(ulist[i], 1)} a.u. L"
	axs[subplot_number].text(1, .99, label, ha='right', va='top', transform=axs[subplot_number].transAxes)
	if subplot_number == 0:
		axs[subplot_number].legend(frameon=False, ncol=1, loc='upper center', bbox_to_anchor=(0, 1, 1, 1))
plt.xlabel("Population Response")
plt.xlim([0,75])
#plt.xticks([-20,0,20],["$\mu$-20","$\mu$","$\mu$+20"])
plt.tight_layout()
plt.savefig(output_path/output_name)
plt.savefig(output_path/(output_name+".svg"), format='svg')
plt.show()
# _____ Figure generation END _____
