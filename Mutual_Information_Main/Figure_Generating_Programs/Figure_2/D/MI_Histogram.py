import numpy as np
import os
from pathlib import Path
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import matplotlib.pyplot as plt

# Generates histogram for mutual information values of single cells for model fit EGFR data

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
# Model mutual information file path
mdl_mi_path = Path("Data/EGFR/Information_Values/mdl_EGFR_nf_1/single_cell_MI.npy")

# Model channel capacity input for ICee
mdl_cc_in_path = Path("Data/EGFR/Information_Values/mdl_EGFR_nf_1/cc_input_distribution.npy")

# Experiment channel capacity input for I
expt_cc_in_path = Path("Data/EGFR/Information_Values/expt_EGFR/cc_input_distribution.npy")

# File path for figure image file
output_path = Path("Figures/Figure_2/D/")

output_file_name = "single_cell_performance"
# _____ File path declarations END _____


# _____ Loading files BEGIN _____
# Loads Mutual information files, information taken from steady state responses
model_single_cell_mi = np.load(mdl_mi_path)

# Loads channel capacity input distributions
mdl_cc_in = np.load(mdl_cc_in_path)
expt_cc_in = np.load(expt_cc_in_path)
# _____ Loading files END _____


# _____ Figure generation BEGIN _____
bin_count = 50
mdl_color = '#0096FF'
expt_pop_color = '#52D452'


# Data range tolerance
drt = 0.025

# Horizontal length of second subplot will be hr times larger than the first
hr = 6

# generates the figure
fig, ax = plt.subplots(1)

# Histogram of model fit single cell performances
ax.hist(model_single_cell_mi,bins=25,color = mdl_color,weights = np.ones_like(model_single_cell_mi)/np.size(model_single_cell_mi),alpha = 0.55)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(axis='y',which = "both",left = False)

# generates inset
axin1 = inset_axes(ax, width='100%', height='100%',bbox_to_anchor=(0.125, .3, .7, .45), bbox_transform=ax.transAxes, loc=2, borderpad=0)

# plots model cc input in inset
axin1.bar(np.arange(np.size(mdl_cc_in))-.2,mdl_cc_in,np.ones(np.size(mdl_cc_in))*.4,alpha = 0.55,color = mdl_color,label = "Model $\mathcal{I}_{Cee}$")
axin1.bar(np.arange(np.size(expt_cc_in))+.2,expt_cc_in,np.ones(np.size(expt_cc_in))*.4,alpha = 0.55,color = expt_pop_color,label = "Expt $\mathcal{I}_{CSA}$")

# removes axis spines
axin1.spines["top"].set_visible(False)
axin1.spines["right"].set_visible(False)
axin1.set_xlabel("CC Input Dose (ng/mL)")

axin1.set_xticks(np.arange(np.size(mdl_cc_in)))
axin1.set_xticklabels(['0','$\\frac{1}{128}$','$\\frac{1}{64}$','$\\frac{1}{32}$','$\\frac{1}{16}$','$\\frac{1}{8}$','$\\frac{1}{4}$','$\\frac{1}{2}$','1','100'])
# sets inset to be in the log scale

ax.set_xlabel("Mutual Information (Bits)")
ax.set_ylabel("Probability")
handles, labels = axin1.get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(.5,1.),ncol = 2,frameon=False)

# generates and saves plot in the tight layout
plt.tight_layout()
plt.savefig(output_path/output_file_name)
plt.savefig(output_path/(output_file_name + ".svg"),format = 'svg')
plt.show()
# _____ Figure generation End _____
