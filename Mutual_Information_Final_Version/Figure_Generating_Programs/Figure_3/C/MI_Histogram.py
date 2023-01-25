import numpy as np
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import Mutual_Information_Final_Version.functions.file_processing_functions.user_input_functions as uif
import matplotlib.pyplot as plt
# Generates histogram for mutual information values of single cells for experiment and model

# sets font for figure text
plt.rcParams.update({'font.size': 15})

# _____ Load arguments BEGIN _____
# Specifies path to file which specifies user specified arguments
user_input_file_path = "/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Figure_Generating_Programs/Figure_3/C/args.txt"
user_arguments = uif.read_file_to_arg_dict(user_input_file_path)
# _____ Load arguments END _____

# _____ File path declarations BEGIN _____
# Experimental mutual information file path
expt_mi_path = user_arguments["expt_mi_path"]

# Model mutual information file path
mdl_mi_path = user_arguments["mdl_mi_path"]

# Experimental channel capacity input file path
expt_cc_in_path = "/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Data/IGFR/Information_Values/Experimental_Step_Dose_MI_CSAR/cc_in.csv"

# Model channel capacity input file path
mdl_cc_in_path = "/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Data/IGFR/Information_Values/Model_CeeMI_Step_Dose/cc_in.csv"

# File path for figure image file
output_path = user_arguments["output_path"]
# _____ File path declarations END _____


# _____ Loading files BEGIN _____
# Loads Mutual information files, information taken from steady state responses
single_cell_mi = np.loadtxt(expt_mi_path, delimiter=",")
model_single_cell_mi = np.loadtxt(mdl_mi_path,delimiter=",")

# Loads channel capacity input distributions
expt_cc_in = np.loadtxt(expt_cc_in_path, delimiter=",")
mdl_cc_in = np.loadtxt(mdl_cc_in_path,delimiter=",")
# _____ Loading files END _____


# _____ Figure generation BEGIN _____
bin_count = user_arguments["bin_count"]
mdl_color = user_arguments["mdl_color"]
expt_color = user_arguments["expt_color"]

# Data range to show
dr = [1.6,2]

# Data range tolerance
drt = 0.025

# Horizontal length of second subplot will be hr times larger than the first
hr = 6

# generates the figure with 2 axis for the break in the x axis
fig, (ax1,ax2) = plt.subplots(1,2,sharey=True,gridspec_kw={'width_ratios': [1, hr]})

# plots the model fit cell performances for information values above the cut-off
ax2.hist(model_single_cell_mi,weights = np.ones_like(model_single_cell_mi)/np.size(model_single_cell_mi),bins = np.linspace(dr[0],dr[1],bin_count-1),alpha = 0.7,range = (dr[0],dr[1]),label = "Model",color = mdl_color)
# plots the model fit cell performances for information values below the cut-off
ax2.hist(single_cell_mi,weights = np.ones_like(single_cell_mi)/np.size(single_cell_mi),bins = np.linspace(dr[0],dr[1],bin_count-1) ,alpha = 0.40,range = (dr[0],dr[1]),label = "Experiment",color = expt_color)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax2.spines["left"].set_visible(False)
ax2.tick_params(axis='y',which = "both",left = False)

# cuts the x lim off between the ranges of the data desired, with some extra tolerance for appearance's sake
ax2.set_xlim([dr[0]-drt,dr[1]])

# gives the same width of the boxes in the right plot, modified by their relative sizes
ltb_width = (np.linspace(dr[0],dr[1],bin_count+1)[1]-np.linspace(dr[0],dr[1],bin_count+1)[0])*hr

# plots the model fit values for information values below the cut-off
ax1.bar(.25,np.count_nonzero(model_single_cell_mi <= dr[0])/np.size(model_single_cell_mi),ltb_width,alpha = 0.7,color = mdl_color)

# plots the experimental values for information values below the cut-off
ax1.bar(.25,np.count_nonzero(single_cell_mi <= dr[0])/np.size(single_cell_mi),ltb_width,alpha = 0.55,color = expt_color)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax1.set_xticks([.25])
ax1.set_xticklabels(["$\leq$1.6"])

# Ensures the range of the second plot is the same as the range of the first
# This makes it easier to have the box lengths be the same for both subplots
ax1.set_xlim([.25-(dr[1]-dr[0]+drt)/2 ,.25+(dr[1]-dr[0]+drt)/2])

# adds the little slash marks at the axis cut
d = .01
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot((1-5.25*d+.01,1+5.25*d+.01), (-d,+d), **kwargs)
kwargs.update(transform=ax2.transAxes)
ax2.plot((-d,+d), (-d,+d), **kwargs)

# generates inset
axin1 = inset_axes(ax2, width='100%', height='100%',bbox_to_anchor=(0.125, .3, .60, .45),
                   bbox_transform=ax2.transAxes, loc=2, borderpad=0)

# plots experimental values in inset
axin1.bar(np.arange(np.size(expt_cc_in))-.15,expt_cc_in,np.ones(np.size(expt_cc_in))*.3,alpha = 0.55,color = expt_color)

# plots model fit values in inset
axin1.bar(np.arange(np.size(mdl_cc_in))+.15,mdl_cc_in,np.ones(np.size(mdl_cc_in))*.3,alpha = 0.7,color = mdl_color)

# removes axis spines
axin1.spines["top"].set_visible(False)
axin1.spines["right"].set_visible(False)
axin1.set_xlabel("CC Input Dose (pM)")

#axin1.set_xlim([dr[0],dr[1]])
#axin1.set_yticks([0,0.01,0.02])

axin1.set_xticks(np.arange(np.size(expt_cc_in)))
axin1.set_xticklabels(["0","17.5","37","125"])

# sets inset to be in the log scale
#axin1.set_yscale('log')

ax2.set_xlabel("Information (Bits)")
ax2.xaxis.set_label_coords(.38, -0.1)
ax1.set_ylabel("Probability")
#plt.yscale('log', nonpositive='clip')
#plt.ylim([0.0007,1])
#plt.yticks([0.001,0.01,0.1])
handles, labels = ax2.get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(.5,1.),ncol = 2,frameon=False)

# generates and saves plot in the tight layout
plt.tight_layout()
plt.savefig(output_path)
plt.show()
# _____ Figure generation End _____
