import numpy as np
import seaborn as sns
import os
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import Mutual_Information_Main.functions.data_processing_functions.sensitivity_analysis as sa
import random
import scipy.stats as st
# Generates scatter plot for MI vs experimental response range

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
# Experimental multi-dose moments file path
expt_md_moment_path = Path("Data/IGFR/Moments/experimental_moments_single_cell/multidose_moments.npy")

# Model multi-dose moments file path
mdl_md_moment_path = Path("Data/IGFR/Information_Values/Model_CeeMI_Step_Dose/mdl_moments_used.npy")

# Experimental mutual information file path
expt_mi_path = Path("Data/IGFR/Information_Values/Experimental_Step_Dose_Single_Cell_Performances/single_cell_MI.npy")

# Model mutual information file path
mdl_mi_path = Path("Data/IGFR/Information_Values/Model_CeeMI_Step_Dose/single_cell_mi_array.npy")

# Location of response range cartoon
rr_cartoon = Path("Figures/Figure_3/left/range.png")

# File path for figure image file
output_path = Path("Figures/Figure_3/left/MI_vs_Response_Range")
# _____ File path declarations END _____


# _____ Loading files BEGIN _____
# Loads Mutual information files, information taken from steady state responses
single_cell_mi = np.load(expt_mi_path)
model_single_cell_mi = np.load(mdl_mi_path)

# Loads moment files
single_cell_moments = np.load(expt_md_moment_path)
model_single_cell_moments = np.load(mdl_md_moment_path)
# _____ Loading files END _____


# _____ Data Processing BEGIN _____
# Generates response ranges by taking the difference in the mean responses
# from inputs of 0 pM to 125 pM at steady state (60 minutes after dose)
response_range = single_cell_moments[:,0]-single_cell_moments[:,3]
#response_range = np.max(single_cell_moments[:,:4],axis = 1)-np.min(single_cell_moments[:,:4],axis = 1)
model_response_range = model_single_cell_moments[:,0]-model_single_cell_moments[:,3]

var_val_ar = 10**np.linspace(-3,3,20)
pc_list = []

print(st.pearsonr(response_range,single_cell_mi))
print(st.pearsonr(model_response_range,model_single_cell_mi))


bin_count = 16
si_list = []
sensitivity_results = sa.bin_y_wrt_x(model_response_range, model_single_cell_mi,
                                     x_samp=response_range, y_samp=single_cell_mi, nbins=16,
                                     bin_center="center_of_bin_mass")
print(sensitivity_results)
var_mdl = sensitivity_results["total_pop_variance"]
var_expt = sensitivity_results["total_samp_variance"]

var_exp_mdl = sensitivity_results["variance_of_expected_values_pop"]
var_exp_expt = sensitivity_results["variance_of_expected_values_samp"]
var_exp_mdl_expt_bin = sensitivity_results["variance_of_expected_values_pop_expt_bin"]

exp_var_mdl = sensitivity_results["expected_value_of_variance_pop"]
exp_var_expt = sensitivity_results["expected_value_of_variance_samp"]

# stores the expected value and locations of bins and the std of the samples within the bins
mdl_bin_center = sensitivity_results["mdl_bin_centers"]
exp_bin_center = sensitivity_results["expt_bin_centers"]
mdl_information_bins = sensitivity_results['model_bin_heights']
exp_information_bins = sensitivity_results['expt_bin_heights']
mdl_information_bins_std_d = sensitivity_results['model_bin_std_d']
exp_information_bins_std_e = sensitivity_results['expt_bin_std_e']
# _____ Data Processing END _____


# _____ Figure generation BEGIN _____
# sets figure colors
mdl_color = (0,33/255,165/255)
expt_color = (200/255,30/255,2/255)
mdl_trend_color = '#39ff14'
expt_trend_color = "cyan"

# initialize figure
fig, ax = plt.subplots(1)

# sets axis labels
ax.set_xlabel("Response Range")
ax.set_ylabel("Mutual Information (Bits)")

# allows kdeplot to have a standard legend label and visual
blue_patch = mpatches.Patch(color=mdl_color, label='Model',alpha = .6)

# Generates kernel density estimation plot for model data
sns.kdeplot(x=model_response_range, y=model_single_cell_mi,levels = [0.01,0.1,0.5,1],color = mdl_color,bw_adjust = 3,fill=True)


# Generates scatter plot for experimental data
val = plt.scatter(response_range,single_cell_mi,label = "Experiment",color = expt_color)
print(st.spearmanr(response_range,single_cell_mi))
print(st.spearmanr(model_response_range,model_single_cell_mi))
# plots model trend line
lin_1 = plt.plot(mdl_bin_center,mdl_information_bins,label = "Model Trend", color = mdl_trend_color,linewidth = 3)

# plots experimental trend line
lin_2 = plt.errorbar(exp_bin_center,exp_information_bins,yerr=exp_information_bins_std_e,label = "Expt Trend", color = expt_trend_color,linewidth = 3)

print(st.pearsonr(mdl_bin_center, mdl_information_bins))
print(st.pearsonr(exp_bin_center, exp_information_bins))

# generates legend
#ax.legend(handles=[blue_patch,lin_1[0],val,lin_2],loc="upper center",frameon=False,ncol = 2, bbox_to_anchor=(.5,1.25))

# removes plot spines
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# tweaks axis properties for visual effect
ax.set_ylim([1,2.02])
ax.set_yticks([1,1.2,1.4,1.6,1.8,2])
ax.set_xticks([100,200,300,400,500])
ax.set_xlim([90,480])

img = plt.imread(rr_cartoon)
# above legend
#newax = fig.add_axes([0.65,0.45,0.3,0.3], anchor='E', zorder=1)
# left of legend
#newax = fig.add_axes([0.3,0.15,0.3,0.3], anchor='E', zorder=1)
# at legend
newax = fig.add_axes([0.55,0.18,0.35,0.35], anchor='E', zorder=1)
newax.imshow(img)
newax.axis('off')

# generates and saves figure in the tight layout
plt.tight_layout()
plt.savefig(output_path)
plt.savefig(f"{output_path}.svg",format = 'svg')
plt.show()
# _____ Figure generation END _____
