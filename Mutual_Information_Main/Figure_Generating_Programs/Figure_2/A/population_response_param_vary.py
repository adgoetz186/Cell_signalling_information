from mpl_toolkits.axes_grid.inset_locator import inset_axes
import pandas as pd
import numpy as np
import sys
import os
from pathlib import Path
import matplotlib.patches as mpatches
import matplotlib.colors as mc
import matplotlib.pyplot as plt

# Generates plot for toy model system
# Shows:
# Mutual information as a function of average initial receptor count I(theta)
# Distribution of single cell mutual information values (pCeeMI(I)) for different amounts of parameter variability
# Average of single cell mutual information values (ICee) for different amounts of parameter variability
# cell state agnostic response based mutual information values (ICSA for different amounts of parameter variability

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
# Cee_MI for variable average initial receptor count, r0
cmi_r0_varies_filename = Path("Data/Toy_Model/Information_Values/framework_comparison/toy_model_system_cmi_r0.csv")
# MI of CSAR for variable average initial receptor count, r0
mi_r0_varies_filename = Path("Data/Toy_Model/Information_Values/framework_comparison/toy_model_system_mi_r0.csv")
# single cell performance under different initial receptor counts, r0
mi_r0_sc_filename = Path("Data/Figure_Generating_Data/Figure_2/B/single_cell_toy_model_system_cmi_single_cell_r0.csv")
# Location to store the main figure
output_folder = Path("Figures/Figure_2/A/")
# Location name of the figure
output_name_suffix = "information_comparison"
# _____ File path declarations END _____

# _____ Loading files BEGIN _____
mi_r0_sc = pd.read_csv(mi_r0_sc_filename, index_col=0, header=0)
cmi_r0 = pd.read_csv(cmi_r0_varies_filename, index_col=0, header=0)
mi_r0 = pd.read_csv(mi_r0_varies_filename, index_col=0, header=0)
# _____ Loading files End _____


# _____ Main code BEGIN _____

# MAIN TEXT FIGURE
# reads parameters from file
input_cv_list = ["-1", "-1/2", "0"]
s_color = "#99ccff"
ceemi_color = "black"
micsar_color = "blue"
sc_color = "blue"
sa = 0.6
sp = [5, 95]

# creates numpy arrays
mi_r0_sc_array = mi_r0_sc.to_numpy()
r0_cv = cmi_r0.columns.to_numpy(dtype=float)
mi_r0_array = mi_r0.to_numpy()
cmi_r0_array = cmi_r0.to_numpy()

# generates figure
fig, ax = plt.subplots(1)
ind_to_plot = [20, 32, 42]
# plots Cee-MI for different receptor distributions
# for i in range(50):
#    plt.hist(cmi_deg_array[:,i],bins=25)
#    plt.show()
# Places points at Cee-MI where histograms go
ax.scatter([r0_cv[20], r0_cv[32], r0_cv[42]],
           [np.average(cmi_r0_array, axis=0)[20], np.average(cmi_r0_array, axis=0)[32],
            np.average(cmi_r0_array, axis=0)[42]], color=ceemi_color)

# generates axis for histograms
ax2 = ax.twiny()
max_r = 40
min_r = 0
ax2.set_xlim([min_r, max_r])
ax.set_xlim([np.min(r0_cv), np.max(r0_cv)])
display_dose_range_log_10 = np.log10(np.max(r0_cv)) - np.log10(np.min(r0_cv))

# sets bases for histograms
adjusted_doses = []
for d_t_h in ind_to_plot:
	adjusted_doses.append(
		(max_r - min_r) * (np.log10(r0_cv[d_t_h]) - np.log10(np.min(r0_cv))) / display_dose_range_log_10)
# plots histograms
for key_ind in range(len(ind_to_plot)):
	ax2.hist(cmi_r0_array[:, ind_to_plot[key_ind]], orientation="horizontal", bottom=adjusted_doses[key_ind], alpha=1,
	         density='pdf', bins=40, color=s_color, zorder=1)
ax2.set_zorder(1)
ax.set_zorder(2)
ax.set_facecolor("none")
lin1 = ax.plot(r0_cv, np.average(cmi_r0_array, axis=0), color=ceemi_color, label="$\mathcal{I}_{Cee}$", zorder=2)
# plots cell spread for different receptor distributions
# ax.fill_between(r0_cv,np.percentile(cmi_deg_array,sp[1],axis=0),np.percentile(cmi_deg_array,sp[0],axis=0),alpha = sa,color = s_color,label = "SC Spread")
# plots MI of CSAR for different receptor distributions
blue_patch = mpatches.Patch(color=s_color, label="$p_{CeeMI}(\mathcal{I})$")
lin2 = ax.plot(r0_cv, np.transpose(mi_r0_array), linestyle="dashed", color=micsar_color, label="$\mathcal{I}_{CSA}$")
ax2.set_xticks([])
# makes x axis in log scale
ax.set_xscale('log')

# sets x,y labels
ax.set_xlabel("Variation of $R_0$ (CV)")
ax.set_ylabel("Mutual Information (Bits)")

# sets y limit between 0 and 4 bits
ax.set_ylim([0, 4])

# removes right and upper axis spines
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)

# Generates inset of single cell performance
axin1 = inset_axes(ax, width='100%', height='100%', bbox_to_anchor=(0.15, .2, .3, .3),
                   bbox_transform=ax.transAxes, loc=2, borderpad=0)

# plots single cell performance in inset
axin1.semilogx(10 ** np.array(mi_r0_sc.columns, dtype=float), -1 * mi_r0_sc_array[-1], label=input_cv_list[-1],
               color=sc_color)

# removes right and upper axis spines in inset
axin1.spines["top"].set_visible(False)
axin1.spines["right"].set_visible(False)

# sets x,y labels for inset
axin1.set_xlabel("$R_0$")
axin1.set_ylabel("$\mathcal{I}(\\theta)$")

# sets y ticks to only be at 0 and 3 bits for inset
axin1.set_yticks([0, 3])

# makes the legend in the upper center of the plot
ax.legend(handles=[lin1[0], lin2[0], blue_patch], loc="upper center", bbox_to_anchor=(.5, 1.1), ncol=3, frameon=False)
# generates and saves figure image in the tight layout
plt.tight_layout()
plt.savefig(output_folder/output_name_suffix)
svg_suff = output_name_suffix+".svg"
plt.savefig(output_folder/svg_suff, format='svg')
plt.show()