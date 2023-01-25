from mpl_toolkits.axes_grid.inset_locator import inset_axes
import pandas as pd
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import Mutual_Information_Final_Version.functions.file_processing_functions.user_input_functions as uif

# Generates plot for toy model system
# Shows:
# Mutual information as a function of average initial receptor count
# Distribution of single cell mutual information values for different amounts of parameter variability
# Average of single cell mutual information values for different amounts of parameter variability
# cell state agnostic response based mutual information values for different amounts of parameter variability

# sets font for figure text
plt.rcParams.update({'font.size': 15})

# _____ Load arguments BEGIN _____
# Specifies path to file which specifies user specified arguments
user_input_file_path = "/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Figure_Generating_Programs/Figure_2/B/args.txt"
user_arguments = uif.read_file_to_arg_dict(user_input_file_path)
# _____ Load arguments END _____

# _____ File path declarations BEGIN _____
# Cee_MI for variable degradation rate, kdeg
cmi_deg_varies_filename = user_arguments["cmi_deg_varies_filename"]
# MI of CSAR for variable degradation rate, kdeg
mi_deg_varies_filename = user_arguments["mi_deg_varies_filename"]
# Cee_MI for variable average initial receptor count, r0
cmi_r0_varies_filename = user_arguments["cmi_r0_varies_filename"]
# MI of CSAR for variable average initial receptor count, r0
mi_r0_varies_filename = user_arguments["mi_r0_varies_filename"]
# single cell performance under different initial receptor counts, r0
mi_r0_sc_filename = user_arguments["mi_r0_sc_filename"]
# single cell performance under different degradation rates, kdeg
mi_deg_sc_filename = user_arguments["mi_deg_sc_filename"]
# Location to store the figure
output_folder = user_arguments["output_folder"]
# Location name of the figure
output_name_suffix = user_arguments["output_name_suffix"]
# _____ File path declarations END _____

# _____ Loading files BEGIN _____
mi_deg_sc = pd.read_csv(mi_deg_sc_filename,index_col=0,header=0)
mi_r0_sc = pd.read_csv(mi_r0_sc_filename,index_col=0,header=0)
cmi_r0 = pd.read_csv(cmi_r0_varies_filename,index_col=0,header=0)
mi_r0 = pd.read_csv(mi_r0_varies_filename,index_col=0,header=0)
# _____ Loading files End _____


# _____ Main code BEGIN _____

# MAIN TEXT FIGURE
# reads parameters from file
input_cv_list = user_arguments["input_cv_list"]
s_color = user_arguments["s_color"]
ceemi_color = user_arguments["ceemi_color"]
micsar_color = user_arguments["micsar_color"]
sc_color = user_arguments["sc_color"]
sa = user_arguments["sa"]
sp = user_arguments["sp"]

# creates numpy arrays
mi_r0_sc_array = mi_r0_sc.to_numpy()
r0_cv = cmi_r0.columns.to_numpy(dtype=float)
mi_r0_array = mi_r0.to_numpy()
cmi_deg_array = cmi_r0.to_numpy()

# generates figure
fig, ax = plt.subplots(1)
ind_to_plot = [20,32,42]
# plots Cee-MI for different receptor distributions
#for i in range(50):
#    plt.hist(cmi_deg_array[:,i],bins=25)
#    plt.show()
# Places points at Cee-MI where histograms go
ax.scatter([r0_cv[20],r0_cv[32],r0_cv[42]],[np.average(cmi_deg_array,axis=0)[20],np.average(cmi_deg_array,axis=0)[32],np.average(cmi_deg_array,axis=0)[42]],color = ceemi_color)

# generates axis for histograms
ax2=ax.twiny()
max_r = 30
min_r = 0
ax2.set_xlim([min_r,max_r])
ax.set_xlim([np.min(r0_cv),np.max(r0_cv)])
display_dose_range_log_10 = np.log10(np.max(r0_cv)) - np.log10(np.min(r0_cv))

# sets bases for histograms
adjusted_doses = []
for d_t_h in ind_to_plot:
	adjusted_doses.append((max_r-min_r)*(np.log10(r0_cv[d_t_h])-np.log10(np.min(r0_cv)))/display_dose_range_log_10)

# plots histograms
for key_ind in range(len(ind_to_plot)):
	ax2.hist(cmi_deg_array[:,ind_to_plot[key_ind]],orientation = "horizontal",bottom = adjusted_doses[key_ind],alpha = 0.5,density='pdf',bins=40,color = s_color)

lin1 = ax.plot(r0_cv,np.average(cmi_deg_array,axis=0),color = ceemi_color,label = "Avg. CeeMI")
# plots cell spread for different receptor distributions
#ax.fill_between(r0_cv,np.percentile(cmi_deg_array,sp[1],axis=0),np.percentile(cmi_deg_array,sp[0],axis=0),alpha = sa,color = s_color,label = "SC Spread")
# plots MI of CSAR for different receptor distributions
blue_patch = mpatches.Patch(color=s_color,label = "p(I)",alpha = .5)
lin2 = ax.plot(r0_cv, np.transpose(mi_r0_array), linestyle = "dashed",color = micsar_color,label = "CSAMI")
ax2.set_xticks([])
# makes x axis in log scale
ax.set_xscale('log')

# sets x,y labels
ax.set_xlabel("Variation of $R_0$")
ax.set_ylabel("Mutual Information (Bits)")

# sets y limit between 0 and 4 bits
ax.set_ylim([0,4])

# removes right and upper axis spines
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)

# Generates inset of single cell performance
axin1 = inset_axes(ax, width='100%', height='100%',bbox_to_anchor=(0.15, .2, .3, .3),
                   bbox_transform=ax.transAxes, loc=2, borderpad=0)

# plots single cell performance in inset
axin1.semilogx(10**np.array(mi_r0_sc.columns,dtype=float),-1*mi_r0_sc_array[-1],label = input_cv_list[-1],color = sc_color)

# removes right and upper axis spines in inset
axin1.spines["top"].set_visible(False)
axin1.spines["right"].set_visible(False)

# sets x,y labels for inset
axin1.set_xlabel("$R_0$")
axin1.set_ylabel("MI (Bits)")

# sets y ticks to only be at 0 and 3 bits for inset
axin1.set_yticks([0,3])

# makes the legend in the upper center of the plot
ax.legend(handles=[lin1[0],lin2[0],blue_patch],loc="upper center", bbox_to_anchor=(.5,1.1),ncol = 3,frameon=False)
# generates and saves figure image in the tight layout
plt.tight_layout()
plt.savefig(f"{output_folder}{output_name_suffix}")
plt.show()

# SUPPLEMENTARY FIGURE

# creates numpy arrays
mi_deg_sc_array = mi_deg_sc.to_numpy()
cmi_deg = pd.read_csv(cmi_deg_varies_filename,index_col=0,header=0)
mi_deg = pd.read_csv(mi_deg_varies_filename,index_col=0,header=0)
mi_deg_array = mi_deg.to_numpy()
cmi_deg_array = cmi_deg.to_numpy()
deg_cv = cmi_deg.columns.to_numpy(dtype=float)

# generates figure
fig, ax = plt.subplots(1)


# Places points at Cee-MI where histograms go
ax.scatter([r0_cv[20],r0_cv[32],r0_cv[42]],[np.average(cmi_deg_array,axis=0)[20],np.average(cmi_deg_array,axis=0)[32],np.average(cmi_deg_array,axis=0)[42]],color = ceemi_color)

# generates axis for histograms
ax2=ax.twiny()
max_r = 1000
min_r = 0
ax2.set_xlim([min_r,max_r])
ax.set_xlim([np.min(r0_cv),np.max(r0_cv)])
display_dose_range_log_10 = np.log10(np.max(r0_cv)) - np.log10(np.min(r0_cv))

# sets bases for histograms
adjusted_doses = []
for d_t_h in ind_to_plot:
	adjusted_doses.append((max_r-min_r)*(np.log10(r0_cv[d_t_h])-np.log10(np.min(r0_cv)))/display_dose_range_log_10)

# plots histograms
for key_ind in range(len(ind_to_plot)):
	ax2.hist(cmi_deg_array[:,ind_to_plot[key_ind]],orientation = "horizontal",bottom = adjusted_doses[key_ind],alpha = 0.5,density='pdf',bins=40,color = s_color)
blue_patch = mpatches.Patch(color=s_color,label = "p(I)",alpha = .5)
ax2.set_xticks([])
# plots Cee-MI for different degradation rate distributions
lin1 = ax.plot(deg_cv,np.average(cmi_deg_array,axis=0),color = ceemi_color,label = "Avg. CeeMI")
# plots cell spread for different degradation rate distributions
#ax.fill_between(deg_cv,np.percentile(cmi_deg_array,sp[1],axis=0),np.percentile(cmi_deg_array,sp[0],axis=0),alpha = sa,color = s_color,label = "SC Spread")
# plots MI of CSAR for different degradation rate distributions
lin2 = ax.plot(deg_cv, np.transpose(mi_deg_array), linestyle = "dashed",color = micsar_color,label = "CSAMI")

# makes x axis in log scale
ax.set_xscale("log")

# sets x,y labels
ax.set_xlabel("Variation of $k_{deg}$")
ax.set_ylabel("Mutual Information (Bits)")

# sets y limit between 0 and 1.75 bits
ax.set_ylim([1.2,1.5])

# removes right and upper axis spines
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)

# Generates inset of single cell performance
axin1 = inset_axes(ax, width='100%', height='100%',bbox_to_anchor=(0.2, .2, .3, .3),
                   bbox_transform=ax.transAxes, loc=2, borderpad=0)

# plots single cell performance in inset
axin1.semilogx(10**np.array(mi_deg_sc.columns,dtype=float),-1*mi_deg_sc_array[-1],label = input_cv_list[-1],color = sc_color)

# removes right and upper axis spines in inset
axin1.spines["top"].set_visible(False)
axin1.spines["right"].set_visible(False)

# sets x,y labels for inset
axin1.set_xlabel("$k_{deg}$")
axin1.set_ylabel("MI (Bits)")

# sets y limits to be at 0 and 1.6 bits for inset
axin1.set_ylim([0,1.6])

# sets y ticks to only be at 0 and 1.5 bits for inset
axin1.set_yticks([0,1.5])

# makes the legend in the upper center of the plot
ax.legend(handles = [lin1[0],lin2[0],blue_patch],loc="upper center", bbox_to_anchor=(.5,1.12),ncol = 3,frameon=False)

# generates and saves figure image in the tight layout
plt.tight_layout()
plt.savefig(f"{output_folder}Supp_{output_name_suffix}")
plt.show()
# _____ Main code END _____
