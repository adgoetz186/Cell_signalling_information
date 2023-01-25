import numpy as np
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import Mutual_Information_Final_Version.functions.file_processing_functions.user_input_functions as uif
import matplotlib.pyplot as plt

# Generates histogram for mutual information values of single cells for model fit EGFR data

# sets font for figure text
plt.rcParams.update({'font.size': 15})

# _____ Load arguments BEGIN _____
# Specifies path to file which specifies user specified arguments
user_input_file_path = "/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Figure_Generating_Programs/Figure_2/F/args.txt"
user_arguments = uif.read_file_to_arg_dict(user_input_file_path)
# _____ Load arguments END _____

# _____ File path declarations BEGIN _____
# Model mutual information file path
mdl_mi_path = user_arguments["mdl_mi_path"]

# Model channel capacity input
mdl_cc_in_path = user_arguments["mdl_cc_in_path"]

# File path for figure image file
output_path = user_arguments["output_path"]
# _____ File path declarations END _____


# _____ Loading files BEGIN _____
# Loads Mutual information files, information taken from steady state responses
model_single_cell_mi = np.loadtxt(mdl_mi_path,delimiter=",")

# Loads channel capacity input distribution
mdl_cc_in = np.loadtxt(mdl_cc_in_path,delimiter=",")
# _____ Loading files END _____


# _____ Figure generation BEGIN _____
bin_count = user_arguments["bin_count"]
mdl_color = user_arguments["mdl_color"]
expt_color = user_arguments["expt_color"]


# Data range tolerance
drt = 0.025

# Horizontal length of second subplot will be hr times larger than the first
hr = 6

# generates the figure
fig, ax = plt.subplots(1)

# Histogram of model fit single cell performances
ax.hist(model_single_cell_mi,bins=25,color = expt_color,weights = np.ones_like(model_single_cell_mi)/np.size(model_single_cell_mi))

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(axis='y',which = "both",left = False)

# generates inset
axin1 = inset_axes(ax, width='100%', height='100%',bbox_to_anchor=(0.125, .3, .7, .45), bbox_transform=ax.transAxes, loc=2, borderpad=0)

# plots model cc input in inset
axin1.bar(np.arange(np.size(mdl_cc_in)),mdl_cc_in,np.ones(np.size(mdl_cc_in))*.8,alpha = 0.7,color = expt_color)

# removes axis spines
axin1.spines["top"].set_visible(False)
axin1.spines["right"].set_visible(False)
axin1.set_xlabel("CC Input Dose (nM)")

axin1.set_xticks(np.arange(np.size(mdl_cc_in)))
axin1.set_xticklabels(['0','$\\frac{1}{128}$','$\\frac{1}{64}$','$\\frac{1}{32}$','$\\frac{1}{16}$','$\\frac{1}{8}$','$\\frac{1}{4}$','$\\frac{1}{2}$','1','100'])
# sets inset to be in the log scale

ax.set_xlabel("Information (Bits)")
ax.set_ylabel("Probability")
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(.5,1.),ncol = 2,frameon=False)

# generates and saves plot in the tight layout
plt.tight_layout()
plt.savefig(output_path)
plt.show()
# _____ Figure generation End _____
