import numpy as np
import seaborn as sns
import random
import os
import scipy.stats as st
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import Mutual_Information_Main.functions.data_processing_functions.sensitivity_analysis as sa
# Generates scatter plot for MI vs initial nuclear foxo

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
# Experimental multi-dose trajectories file path
expt_md_moment_path = Path("Data/IGFR/experimental_data/adjusted_multidose/multi_dose.npy")

# Model multi-dose moments file path
mdl_md_moment_path = Path("Data/IGFR/Information_Values/Model_CeeMI_Step_Dose/mdl_moments_used.npy")

# Experimental mutual information file path
expt_mi_path = Path("Data/IGFR/Information_Values/Experimental_Step_Dose_Single_Cell_Performances/single_cell_MI.npy")

# Model mutual information file path
mdl_mi_path = Path("Data/IGFR/Information_Values/Model_CeeMI_Step_Dose/single_cell_mi_array.npy")

# Location of initial Foxo cartoon
ifo_cartoon = Path("Figures/Figure_3/right/Initial_n_Foxo_Cartoon.png")

# File path for figure image file
output_path = Path("Figures/Figure_3/right/")

# Name of output
output_name = "MI_vs_IFO"
# _____ File path declarations END _____


# _____ Loading files BEGIN _____
# Loads Mutual information files, information taken from steady state responses
single_cell_mi = np.load(expt_mi_path)
model_single_cell_mi = np.load(mdl_mi_path)

# Loads experimental trajectory file
single_cell_nuclear_foxo = np.load(expt_md_moment_path)

# Loads model moment file
model_single_cell_moments = np.load(mdl_md_moment_path)
# _____ Loading files END _____


# _____ Data Processing BEGIN _____
# Finds each cells initial nuclear foxo averaged over all recorded pre dose measurements for each cell
single_cell_mean_initial_nuclear_foxo = np.average(single_cell_nuclear_foxo[1:,:11],axis = 1)

# Model initial nuclear foxo
mdl_ifoxo = model_single_cell_moments[:,0]

bin_count = 16
sensitivity_results = sa.bin_y_wrt_x(mdl_ifoxo,model_single_cell_mi,x_samp = single_cell_mean_initial_nuclear_foxo,y_samp = single_cell_mi, nbins = 16, bin_center = "center_of_bin_mass")

var_mdl = sensitivity_results["total_pop_variance"]
var_expt = sensitivity_results["total_samp_variance"]

var_exp_mdl = sensitivity_results["variance_of_expected_values_pop"]
var_exp_expt = sensitivity_results["variance_of_expected_values_samp"]

exp_var_mdl = sensitivity_results["expected_value_of_variance_pop"]
exp_var_expt = sensitivity_results["expected_value_of_variance_samp"]


# Finds bin edges while ensuring each bin gets a similar number of points
percentile_edges = np.linspace(0, 100, bin_count + 1)
mdl_bins = np.percentile(mdl_ifoxo, percentile_edges)
exp_bins = np.percentile(single_cell_mean_initial_nuclear_foxo, percentile_edges)

# Finds center of bins for plotting
mdl_bin_center = (mdl_bins[1:] + mdl_bins[:-1]) / 2
exp_bin_center = (exp_bins[1:] + exp_bins[:-1]) / 2

# This adds a little to the upper bounds of the far right bin to avoid having a single point in the highest bin
mdl_bins[-1] += 1
exp_bins[-1] += 1

# Finds which bin each point belongs to, binning is done based on initial nuclear foxo
mdl_ifoxo_digi = np.digitize(mdl_ifoxo, mdl_bins)
exp_ifoxo_digi = np.digitize(single_cell_mean_initial_nuclear_foxo, exp_bins)

# For the standard error bins the model using the same bins as used by the experimental data for comparison
mdl_ifoxo_exp_bin_digi = np.digitize(mdl_ifoxo, exp_bins)

# initializes numpy arrays to store bin data
mdl_information_bins = np.zeros(bin_count)
exp_information_bins = np.zeros(bin_count)
mdl_information_bins_std = np.zeros(bin_count)
exp_information_bins_std = np.zeros(bin_count)

for i in range(1, bin_count + 1):
	# Calculates the binned average information values
	mdl_information_bins[i - 1] = np.average(model_single_cell_mi[np.argwhere(mdl_ifoxo_digi == i)])
	exp_information_bins[i - 1] = np.average(single_cell_mi[np.argwhere(exp_ifoxo_digi == i)])
	
	# Calculates the standard error for the sample data
	mdl_information_bins_std[i - 1] = np.std(model_single_cell_mi[np.argwhere(mdl_ifoxo_exp_bin_digi == i)])
	exp_information_bins_std[i - 1] = mdl_information_bins_std[i - 1] / np.size(
		single_cell_mi[np.argwhere(exp_ifoxo_digi == i)]) ** .5
# _____ Data Processing END _____


# _____ Figure generation BEGIN _____
# sets the figure colors
mdl_color = (0,33/255,165/255)
expt_color = (200/255,30/255,2/255)
mdl_trend_color = '#39ff14'
expt_trend_color = "cyan"


# initialize figure
fig, ax = plt.subplots(1)

# sets axis labels
ax.set_xlabel("Initial nuclear Foxo")
ax.set_ylabel("Mutual Information (Bits)")

# allows kdeplot to have a standard legend label and visual
blue_patch = mpatches.Patch(color=mdl_color, label='Model', alpha=.6)

# Generates kernel density estimation plot for model data
sns.kdeplot(x=mdl_ifoxo, y=model_single_cell_mi, levels=[0.01, 0.1, 0.5, 1], color=mdl_color,bw_adjust=3, fill=True, alpha=1)



# Generates scatter plot for experimental data
val_2 = plt.scatter(single_cell_mean_initial_nuclear_foxo,single_cell_mi,label = "Experiment",color = expt_color)

# removes plot spines
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# plots model trend line
lin_1 = plt.plot(mdl_bin_center,mdl_information_bins,label = "Model Trend", color = mdl_trend_color,linewidth = 3)

# plots experimental trend line
lin_2 = plt.errorbar(exp_bin_center,exp_information_bins,yerr=exp_information_bins_std,label = "Expt Trend", color = expt_trend_color,linewidth = 3)

# generates legend

img = plt.imread(ifo_cartoon)

newax = fig.add_axes([0.55,0.18,0.35,0.35], anchor='E', zorder=1)
newax.imshow(img)
newax.axis('off')
# tweaks axis properties for visual effect
ax.set_ylim([1,2.02])
ax.set_xlim([280,780])

plt.tight_layout()
plt.savefig(output_path / output_name)
plt.savefig(output_path / (output_name + ".svg"),format='svg')
plt.show()
fig,ax = plt.subplots(1)
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.spines["bottom"].set_visible(False)
ax.legend(handles=[blue_patch,lin_1[0],val_2,lin_2], loc="center",frameon=False,ncol=2)
plt.savefig(output_path / (output_name + "_legend"))
plt.savefig(output_path / (output_name + "_legend.svg"),format='svg')
plt.show()
# _____ Figure generation END _____
