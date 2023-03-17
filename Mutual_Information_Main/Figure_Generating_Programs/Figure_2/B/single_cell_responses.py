import numpy as np
import time
import sys
import os
from pathlib import Path
import scipy.stats as st
import matplotlib.pyplot as plt
from Mutual_Information_Main.functions.toy_model_functions import Toy_Model_Functions as tmf
import pandas as pd

# Generates schematic showing cell state agnostic responses and single cell responses

plt.rcParams.update({'font.size': 15})
plt.rcParams["hatch.linewidth"] = 4

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
# location of experimental CSAR MI
output_dir = Path("Figures/Figure_2/B")
output_name_suffix = "response_distributions"
# _____ File path declarations END _____

# _____ Figure generation BEGIN _____
# Coefficient of variation of input distribution
input_cv = 1

# Mean of input distribution
input_mean = 10

# Number of points to plot for inset
parameter_partitions = 1000

# average degradation rate, k_deg
deg_mean = 5

# average initial receptor count value, r0
r0_mean = 500

# rate of binding, k_bind
k_bind = 1

# rate of unbinding, k_unbind
k_unbind = 10

# max allowed response, responses above will be truncated. Should be significantly higher than r0 mean.
max_response = 25000

# input signal percentiles used to generate responses for R0 varying population
in_p_r0_var = np.array([45,55])/100

# input signal percentiles used to generate responses for kdeg varying population
in_p_deg_var = np.array([35,65])/100

# parameter percentiles used to generate single cell responses
param_p = np.array([1,99])/100

# cv of r0 parameter to use for case where r0 varies
r0_cv = 10**-.75

# creates dataframes to store results
cmir0 = pd.DataFrame(columns=["cell_mi", "input_cv", "parameter_cv"])
cmideg = pd.DataFrame(columns=["cell_mi", "input_cv", "parameter_cv"])

# declares colors for plot
low_in_color = (0/255,33/255,165/255,.2)
high_in_color = (250 / 255, 70 / 255, 22 / 255, 0.2)
overlap_colors = [(0 / 255, 33 / 255, 165 / 255, .5),(250 / 255, 70 / 255, 22 / 255, 1)]

# sets plot limits for R0 varying population
xlim_r0 = [0,400]
ylim_r0 = [0,0.045]

# sets plot limits for kdeg varying population
xlim_deg = [0,50]
ylim_deg = [0,0.1]

# obtains the variance of the input distribution
input_variance = (input_cv * input_mean) ** 2
for uvarind in range(np.size(input_variance)):
	
	# obtains inputs
	u_scale1 = input_variance / input_mean
	uvar_shape1 = input_mean / u_scale1
	ulist = [st.gamma.ppf(in_p_r0_var[0], uvar_shape1, scale=u_scale1),st.gamma.ppf(in_p_r0_var[1], uvar_shape1, scale=u_scale1)]
	
	# obtains the variance of the r0 distribution
	r0_variance = (r0_mean * r0_cv) ** 2
	
	# obtains the binned values of r0 to use
	r0_scale1 = r0_variance / r0_mean
	r0var_shape1 = r0_mean / r0_scale1
	r0edge = np.linspace(0, 1, parameter_partitions + 1)
	r0center = (r0edge[1:] + r0edge[:-1]) / 2
	r0list = st.gamma.ppf(r0center, r0var_shape1, scale=r0_scale1)
	
	# obtains the distribution of responses for both single cell systems
	cond_resp_1 = tmf.pxgut(ulist, deg_mean, st.gamma.ppf(param_p[0], r0var_shape1, scale=r0_scale1), k_bind, k_unbind, 5)
	cond_resp_2 = tmf.pxgut(ulist, deg_mean, st.gamma.ppf(param_p[1], r0var_shape1, scale=r0_scale1), k_bind, k_unbind, 5)
	
	# obtains the distribution of responses for the CSAR
	cond_resp = np.zeros_like(tmf.pxgut(ulist, deg_mean, r0list[-1], k_bind, k_unbind, 5))
	for r0_val in r0list:
		cond_resp[:,:(1+int(np.floor(r0_val*5)))] += tmf.pxgut(ulist, deg_mean, r0_val, k_bind, k_unbind, 5)
	for i in range(np.shape(cond_resp)[0]):
		cond_resp[i] /= np.sum(cond_resp[i])
		
	# generates plot
	fig, (ax1,ax2) = plt.subplots(2,1,sharex = True)
	fig.supylabel("Probability")
	
	# plots CSARs and overlap
	ax1.fill(np.arange(0, np.shape(cond_resp)[1]), cond_resp[0],fc = low_in_color,ec = "black", label = "Low Signal")
	ax1.fill(np.arange(0, np.shape(cond_resp)[1]), cond_resp[1],fc = high_in_color,ec = "black", label = "High Signal")
	ax1.fill(np.arange(0, np.shape(cond_resp)[1]), np.minimum(cond_resp[0], cond_resp[1]) ,fc=(1, 1, 1, 1))
	ax1.fill(np.arange(0, np.shape(cond_resp)[1]), np.minimum(cond_resp[0],cond_resp[1]) , fc= overlap_colors[0],hatch="//",ec=overlap_colors[1],lw=0, label = "Overlap")
	
	# sets limits
	ax1.set_ylim([ylim_r0[0], ylim_r0[1]])
	ax1.set_xlim([xlim_r0[0], xlim_r0[1]])
	
	# removes axis spines
	ax1.spines["top"].set_visible(False)
	ax1.spines["right"].set_visible(False)
	
	# draws legend
	ax1.legend(frameon=False)
	
	# adds title
	ax1.set_title("CSAR")

	# plots cell A responses and overlap
	ax2.fill(np.arange(0, np.shape(cond_resp_1)[1]), cond_resp_1[0] , fc=low_in_color, ec="black")
	ax2.fill(np.arange(0, np.shape(cond_resp_1)[1]), cond_resp_1[1] , fc=high_in_color, ec="black")
	ax2.fill(np.arange(0, np.shape(cond_resp_1)[1]), np.minimum(cond_resp_1[0], cond_resp_1[1]) ,fc=(1, 1, 1, 1))
	ax2.fill(np.arange(0, np.shape(cond_resp_1)[1]), np.minimum(cond_resp_1[0], cond_resp_1[1]) ,fc=overlap_colors[0], hatch="//",ec=overlap_colors[1], lw=0)
	
	# plots cell B responses and overlap
	ax2.fill(np.arange(0, np.shape(cond_resp_2)[1]), cond_resp_2[0] , fc=low_in_color, ec="black")
	ax2.fill(np.arange(0, np.shape(cond_resp_2)[1]), cond_resp_2[1] , fc=high_in_color, ec="black")
	ax2.fill(np.arange(0, np.shape(cond_resp_2)[1]), np.minimum(cond_resp_2[0], cond_resp_2[1]), fc=(1, 1, 1, 1))
	ax2.fill(np.arange(0, np.shape(cond_resp_2)[1]), np.minimum(cond_resp_2[0], cond_resp_2[1]), fc=overlap_colors[0], hatch="//", ec=overlap_colors[1], lw=0)
	
	
	ax2.set_xlabel("Response")
	
	# adds title
	ax2.set_title("Single Cell")
	
	# sets plot limits
	ax2.set_ylim([ylim_r0[0],ylim_r0[1]])
	
	# removes spines
	ax2.spines["top"].set_visible(False)
	ax2.spines["right"].set_visible(False)
	
	# generates and saves plot in tight layout
	plt.tight_layout()
	plt.savefig(output_dir/output_name_suffix)
	plt.savefig(output_dir/(output_name_suffix + ".svg"),format='svg')
	plt.show()
	# The following section follows the same setup as above, only kdeg and r0 are switched
	ulist = [st.gamma.ppf(in_p_deg_var[0], uvar_shape1, scale=u_scale1), st.gamma.ppf(in_p_deg_var[1], uvar_shape1, scale=u_scale1)]
	r0_val = 50
	k_deg_cv = 1
	deg_variance = (deg_mean * k_deg_cv) ** 2
	deg_scale1 = deg_variance / deg_mean
	degvar_shape1 = deg_mean / deg_scale1
	degedge = np.linspace(0, 1, parameter_partitions + 1)
	degcenter = (degedge[1:] + degedge[:-1]) / 2
	deglist = st.gamma.ppf(degcenter, degvar_shape1, scale=deg_scale1)
	cond_resp_1 = tmf.pxgut(ulist,st.gamma.ppf(param_p[0], degvar_shape1, scale=deg_scale1),r0_val,k_bind,k_unbind,5)
	cond_resp_2 = tmf.pxgut(ulist, st.gamma.ppf(param_p[1], degvar_shape1, scale=deg_scale1), r0_val, k_bind, k_unbind, 5)
	cond_resp = np.zeros_like(cond_resp_1)
	for deg_val in deglist:
		cond_resp += tmf.pxgut(ulist,deg_val,r0_val,k_bind,k_unbind,5)
	for i  in range(np.shape(cond_resp)[0]):
		cond_resp[i] /= np.sum(cond_resp[i])
	
	fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
	fig.supylabel("Probability")
	print(cond_resp[0])
	print(cond_resp[1])
	ax1.fill_between(np.arange(0, np.shape(cond_resp)[1]), cond_resp[0],np.zeros_like(cond_resp[0]), fc=low_in_color, ec="black",label="Low Signal")
	ax1.fill_between(np.arange(0, np.shape(cond_resp)[1]), cond_resp[1],np.zeros_like(cond_resp[0]), fc=high_in_color, ec="black",label="High Signal")
	ax1.fill_between(np.arange(0, np.shape(cond_resp)[1]), np.minimum(cond_resp[0], cond_resp[1]),np.zeros_like(cond_resp[0]),fc=(1, 1, 1, 1))
	ax1.fill_between(np.arange(0, np.shape(cond_resp)[1]), np.minimum(cond_resp[0], cond_resp[1]),np.zeros_like(cond_resp[0]),fc=overlap_colors[0], hatch="//",ec=overlap_colors[1], lw=0, label="Overlap")
	#ax1.set_ylim([ylim[0], ylim[1]])
	ax1.set_xlim([xlim_deg[0], xlim_deg[1]])
	ax1.spines["top"].set_visible(False)
	ax1.spines["right"].set_visible(False)
	ax1.legend(frameon=False)
	ax1.set_title("CSAR")
	
	ax2.fill_between(np.arange(0, np.shape(cond_resp_1)[1]), cond_resp_1[0],np.zeros_like(cond_resp_1[0]), fc=low_in_color, ec="black")
	ax2.fill_between(np.arange(0, np.shape(cond_resp_1)[1]), cond_resp_1[1],np.zeros_like(cond_resp_1[1]), fc=high_in_color, ec="black")
	ax2.fill_between(np.arange(0, np.shape(cond_resp_1)[1]), np.minimum(cond_resp_1[0], cond_resp_1[1]),np.zeros_like(cond_resp_1[1]),fc=(1, 1, 1, 1))
	ax2.fill_between(np.arange(0, np.shape(cond_resp_1)[1]), np.minimum(cond_resp_1[0], cond_resp_1[1]),np.zeros_like(cond_resp_1[1]),fc=overlap_colors[0], hatch="//",ec=overlap_colors[1], lw=0)
	
	ax2.fill_between(np.arange(0, np.shape(cond_resp_2)[1]), cond_resp_2[0],np.zeros_like(cond_resp_2[0]), fc=low_in_color,ec="black")
	ax2.fill_between(np.arange(0, np.shape(cond_resp_2)[1]), cond_resp_2[1],np.zeros_like(cond_resp_2[1]), fc=high_in_color,ec="black")
	ax2.fill_between(np.arange(0, np.shape(cond_resp_2)[1]), np.minimum(cond_resp_2[0], cond_resp_2[1]),np.zeros_like(cond_resp_2[1]),fc=(1, 1, 1, 1))
	ax2.fill_between(np.arange(0, np.shape(cond_resp_2)[1]), np.minimum(cond_resp_2[0], cond_resp_2[1]),np.zeros_like(cond_resp_2[1]),fc=overlap_colors[0], hatch="//",ec=overlap_colors[1], lw=0)
	ax2.set_xlabel("Response")
	ax2.set_title("Single Cell")
	#ax2.set_ylim([ylim[0], ylim[1]])
	ax2.spines["top"].set_visible(False)
	ax2.spines["right"].set_visible(False)
	plt.tight_layout()
	#plt.show()
	