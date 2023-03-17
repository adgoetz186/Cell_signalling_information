import numpy as np
import seaborn as sns
import random
import pandas as pd
import scipy.stats as st
import matplotlib.ticker as tk
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import Mutual_Information_Main.functions.data_processing_functions.sensitivity_analysis as sa
from pathlib import Path
import os
plt.rcParams.update({'font.size': 15})

# Generates scatter plot for MI vs initial nuclear foxo

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
# Experimental multi-dose trajectories file path
expt_md_moment_path = "Data/IGFR/experimental_data/adjusted_multidose/multi_dose.npy"

# Path to single cell parameters
params = "Data/IGFR/Information_Values/Model_CeeMI_Step_Dose/mdl_params_used.npy"
params = 10**np.load(params)
# Model multi-dose moments file path
mdl_md_moment_path = "Data/IGFR/Information_Values/Model_CeeMI_Step_Dose/mdl_moments_used.npy"

# Model mutual information file path
mdl_mi_path = "Data/IGFR/Information_Values/Model_CeeMI_Step_Dose/single_cell_mi_array.npy"

# File path for supp figure image file
supp_output_path = "Figures/Supplementary/Additional_Params_Density_Plots"
# _____ File path declarations END _____


# _____ Loading files BEGIN _____
# Loads Mutual information files, information taken from steady state responses
model_single_cell_mi = np.load(mdl_mi_path)

# Loads experimental trajectory file
single_cell_nuclear_foxo = np.load(expt_md_moment_path)

# Loads model moment file
model_single_cell_moments = np.load(mdl_md_moment_path)
# _____ Loading files END _____


# _____ Data Processing BEGIN _____
# Finds each cells initial nuclear foxo averaged over all recorded pre dose measurements for each cell
single_cell_mean_initial_nuclear_foxo = np.average(single_cell_nuclear_foxo[1:,:11],axis = 1)
#single_cell_mean_initial_nuclear_foxo = single_cell_nuclear_foxo[1:,0]
# Model initial nuclear foxo
mdl_ifoxo = model_single_cell_moments[:,0]


bin_count = 16
param_dict = {}
param_dict["$K_{D}$"] = np.log10(params[:,3]/params[:,2])
param_dict["$K_{Rp}$"] = np.log10(params[:,4]/params[:,5])
param_dict["$K_{Ap}$"] = np.log10(params[:,7]/params[:,6])
param_dict["$K_{Fp}$"] = np.log10(params[:,8]/params[:,9])
param_dict["$K_{NT}$"] = np.log10(params[:,11]/params[:,10])
param_dict["Average Initial IGFR"] = np.log10(params[:,0])
param_dict["Total Akt"] = np.log10(params[:,-2])
param_dict["Total Foxo"] = np.log10(params[:,-1])
print(np.average(params[:,-1]))
#dnl_param_list = ["Total Foxo"]
dnl_param_list = []
s_dict = {}
for i in param_dict.keys():
	# the sensitivity analysis algorithm is largely unaffected by log scale, this is only needed to get bin edges
	# which look like in matplotlib's log scale
	if i not in dnl_param_list:
		sar = sa.bin_y_wrt_x(10**param_dict[i],model_single_cell_mi,nbins = 16, bin_center = "center_of_bin_mass")
	else:
		sar = sa.bin_y_wrt_x(param_dict[i],model_single_cell_mi,nbins = 16, bin_center = "center_of_bin_mass")
	s_dict[i] = sar
# List of parameters which are not in log scale

# _____ Data Processing END _____


# _____ Figure generation BEGIN _____
# sets the figure colors
mdl_color = (0,33/255,165/255)
expt_color = (200/255,30/255,2/255)
mdl_trend_color = '#39ff14'
expt_trend_color = "cyan"

#plt.scatter(param_dict["Total Foxo"],param_dict["$K_{pFoxo}$"])
#plt.show()
#plt.scatter(np.log10(param_dict["Total Foxo"]),model_single_cell_mi)
#print(st.pearsonr(np.log10(param_dict["Total Foxo"]), model_single_cell_mi))
#plt.show()
#plt.scatter(param_dict["Total Foxo"],param_dict["Dissociation Rate"])

pearson_results = np.zeros((8,2))
count = 0
fig,axs = plt.subplots(4,2,sharey=True,figsize=(10, 11))
# allows kdeplot to have a standard legend label and visual
blue_patch = mpatches.Patch(color=mdl_color, label='Model', alpha=.6)
for i in param_dict.keys():
	#plt.errorbar(10**sensitivity_results_fio["mdl_bin_centers"],sensitivity_results_fio["model_bin_heights"],yerr=sensitivity_results_fio["model_bin_std_d"],color = mdl_trend_color,linewidth = 3)
	
	if i not in dnl_param_list:
		print("start")
		sns.kdeplot(x=10**param_dict[i], y=model_single_cell_mi, levels=[0.01, 0.1, 0.5, 1], color=mdl_color, bw_adjust=3,fill=True, alpha=1, log_scale=(10, False),ax=axs[count%4,count//4])
		mdl_trend = axs[count%4,count//4].plot(s_dict[i]["mdl_bin_centers"],s_dict[i]["model_bin_heights"],color = mdl_trend_color,linewidth = 3,label="Model Trend")
		min_x_log = np.log10(np.min(s_dict[i]["mdl_bin_centers"]))
		max_x_log = np.log10(np.max(s_dict[i]["mdl_bin_centers"]))
		
		diff = max_x_log-min_x_log
		adj_min = 10**(min_x_log-0.10*diff)
		adj_max = 10 ** (max_x_log + 0.1 * diff)
		axs[count%4,count//4].set_xlim([adj_min, adj_max])
		axs[count%4,count//4].set_xscale("log")
	else:
		axs[count%4,count//4].plot(s_dict[i]["mdl_bin_centers"],s_dict[i]["model_bin_heights"],color = mdl_trend_color,linewidth = 3)
		min_x_log = np.min(s_dict[i]["mdl_bin_centers"])
		max_x_log = np.max(s_dict[i]["mdl_bin_centers"])
		diff = max_x_log-min_x_log
		adj_min = min_x_log-0.05*diff
		adj_max = max_x_log + 0.05 * diff
		axs[count%4,count//4].set_xlim([adj_min, adj_max])
	if i == "Total Foxo":
		axs[count % 4, count // 4].set_xticks([])
	axs[count%4,count//4].set_ylim([1,2.01])
	axs[count%4,count//4].set_xlabel(i)
	pearsonr = st.pearsonr(param_dict[i], model_single_cell_mi)
	pearson_results[count,0] = pearsonr[0]
	pearson_results[count, 1] = pearsonr[1]
	#ax.set_title(f"$\\rho$ = {np.round(pearsonr[0],4)}, log(pval) = {np.round(np.log10(pearsonr[1]),2)}")
	axs[count%4,count//4].spines["top"].set_visible(False)
	axs[count%4,count//4].spines["right"].set_visible(False)
	count += 1
fig.supylabel("Mutual Information (Bits)")
plt.tight_layout()

plt.savefig(f"{supp_output_path}/Scatter_Plots")
plt.savefig(f"{supp_output_path}/Scatter_Plots.svg",format = 'svg')
plt.show()
plt.clf()
fig,ax = plt.subplots(1)
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.spines["bottom"].set_visible(False)
ax.legend(handles=[blue_patch,mdl_trend[0]], loc="center",frameon=False,ncol=2)
plt.savefig(f"{supp_output_path}/Scatter_Plots_Legend.svg",format='svg')
plt.show()


# Generates the table of pearson values, the generated table is then used to create the table displayed in the paper
fig, ax = plt.subplots()
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.spines["bottom"].set_visible(False)

table_dict = {}
table_dict["Parameter"] = list(param_dict.keys())
table_dict["Pearson Cor"] = np.round(pearson_results,3)[:,0]
table_dict["log p value"] = np.log10(pearson_results[:,1])
Table = pd.DataFrame(table_dict)

ax.table(cellText=Table.values,colLabels=Table.columns,loc='center')
plt.savefig(f"{supp_output_path}/Pearson_r.svg")
plt.show()
# _____ Figure generation END _____
