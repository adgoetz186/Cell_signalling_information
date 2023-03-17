import numpy as np
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import os
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
# Generates histogram for mutual information values of single cells for experiment and model

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
# Experimental sc mutual information file path
expt_mi_path = Path("Data/IGFR/Information_Values/Experimental_Step_Dose_Single_Cell_Performances/single_cell_MI.npy")

# Model mutual information file path
mdl_mi_path = Path("Data/IGFR/Information_Values/Model_CeeMI_Step_Dose/single_cell_mi_array.npy")

# Experimental Icee channel capacity input file path
expt_cc_in_path = Path("Data/IGFR/Information_Values/Experimental_Step_Dose_Single_Cell_Performances/cc_input_distribution.npy")

# Experimental I channel capacity input file path
expt_I_cc_in_path = Path("Data/IGFR/Information_Values/Experimental_Step_Dose_Population_Level_Performance/cc_input_distribution.npy")

# Model channel capacity input file path
mdl_cc_in_path = Path("Data/IGFR/Information_Values/Model_CeeMI_Step_Dose/cc_in.npy")

# File path for figure image file
output_path = Path("Figures/Figure_2/E/single_cell_performance")
# _____ File path declarations END _____


# _____ Loading files BEGIN _____
# Loads Mutual information files, information taken from steady state responses
single_cell_mi = np.load(expt_mi_path)
model_single_cell_mi = np.load(mdl_mi_path)

# Loads channel capacity input distributions
expt_cc_in = np.load(expt_cc_in_path)
mdl_cc_in = np.load(mdl_cc_in_path)
expt_I_cc_in = np.load(expt_I_cc_in_path)
# _____ Loading files END _____


# _____ Figure generation BEGIN _____
bin_count = 50
mdl_color = '#0096FF'
expt_color = '#ff85a2'
expt_pop_color = '#52D452'

# Data range to show
dr = [1.6,2]

# Data range tolerance
drt = 0.025

# Horizontal length of second subplot will be hr times larger than the first
hr = 6

# generates the figure with 2 axis for the break in the x axis
fig, (ax1,ax2) = plt.subplots(1,2,sharey=True,gridspec_kw={'width_ratios': [1, hr]})

# plots the model fit cell performances for information values above the cut-off
mdl_sc = ax2.hist(model_single_cell_mi,weights = np.ones_like(model_single_cell_mi)/np.size(model_single_cell_mi),bins = np.linspace(dr[0],dr[1],bin_count-1),alpha = 0.55,range = (dr[0],dr[1]),label = "Model $\mathcal{I}_{Cee}$",color = mdl_color)
# plots the model fit cell performances for information values below the cut-off
expt_sc = ax2.hist(single_cell_mi,weights = np.ones_like(single_cell_mi)/np.size(single_cell_mi),bins = np.linspace(dr[0],dr[1],bin_count-1) ,alpha = 0.55,range = (dr[0],dr[1]),label = 'Expt $\mathcal{I}_{Cee}$',color = expt_color)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax2.spines["left"].set_visible(False)
ax2.tick_params(axis='y',which = "both",left = False)

# cuts the x lim off between the ranges of the data desired, with some extra tolerance for appearance's sake
ax2.set_xlim([dr[0]-drt,dr[1]])

# gives the same width of the boxes in the right plot, modified by their relative sizes
ltb_width = (np.linspace(dr[0],dr[1],bin_count+1)[1]-np.linspace(dr[0],dr[1],bin_count+1)[0])*hr

# plots the model fit values for information values below the cut-off
ax1.bar(.25,np.count_nonzero(model_single_cell_mi <= dr[0])/np.size(model_single_cell_mi),ltb_width,alpha = 0.55,color = mdl_color)

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

# plots experimental sc values in inset
axin1.bar(np.arange(np.size(expt_cc_in)),expt_cc_in,np.ones(np.size(expt_cc_in))*.3,alpha = 0.55,color = expt_color)

# plots experimental pop values in inset
axin1.bar(np.arange(np.size(expt_I_cc_in))+.3,expt_I_cc_in,np.ones(np.size(expt_I_cc_in))*.3,alpha = 0.55,color = expt_pop_color)
expt_pop = mpatches.Patch(color=expt_pop_color, label='Expt $\mathcal{I}$',alpha = .6)


# plots model fit values in inset
axin1.bar(np.arange(np.size(mdl_cc_in))-.3,mdl_cc_in,np.ones(np.size(mdl_cc_in))*.3,alpha = 0.55,color = mdl_color)

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

ax2.set_xlabel("Mutual Information (Bits)")
ax2.xaxis.set_label_coords(.38, -0.1)
ax1.set_ylabel("Probability")
#plt.yscale('log', nonpositive='clip')
#plt.ylim([0.0007,1])
#plt.yticks([0.001,0.01,0.1])
handles, labels = ax2.get_legend_handles_labels()
handles.append(expt_pop)
labels.append('Expt $\mathcal{I}_{CSA}$')
fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(.5,1.),ncol = 3,frameon=False)

# generates and saves plot in the tight layout
plt.tight_layout()
plt.savefig(output_path)
plt.savefig(f"{output_path}.svg",format='svg')
plt.show()
# _____ Figure generation End _____
