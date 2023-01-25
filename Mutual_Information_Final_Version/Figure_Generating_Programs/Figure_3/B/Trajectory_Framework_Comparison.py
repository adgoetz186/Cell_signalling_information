import numpy as np
from matplotlib import pyplot as plt
import Mutual_Information_Final_Version.functions.file_processing_functions.user_input_functions as uif

# Generates plot for comparing cell state agnostic channel with cell state dependent channel for IGF-Akt-FOXO system
# Also provides channel capacity input distributions

# sets font for figure text
plt.rcParams.update({'font.size': 15})

# _____ Load arguments BEGIN _____
# Specifies path to file which specifies user specified arguments
user_input_file_path = "/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Figure_Generating_Programs/Figure_3/B/args.txt"
user_arguments = uif.read_file_to_arg_dict(user_input_file_path)
# _____ Load arguments END _____


# _____ File path declarations BEGIN _____
# CeeMI array pathway
ceemi_path = user_arguments["ceemi_path"]

# CeeMI cc input pathway
ceemi_input_path = user_arguments["ceemi_input_path"]

# expt MI pathway
expt_mi_path = user_arguments["expt_mi_path"]

# expt MI input pathway
expt_mi_input_path = user_arguments["expt_mi_input_path"]

# mi for single cells pathway
mi_single_cell_path = user_arguments["mi_single_cell_path"]

# File path for figure image file
output_path = user_arguments["output_path"]
# _____ File path declarations END _____


# _____ Loading files BEGIN _____
# Loads cell state conditioned mutual information trajectories
ceemi = np.loadtxt(ceemi_path,delimiter=",")

# Loads cell state conditioned mutual information input at channel capacity
cc_ceemi_input = np.loadtxt(ceemi_input_path,delimiter=",")

# Loads experimental mutual information trajectories
exp_mi = np.loadtxt(expt_mi_path,delimiter=",")

# Loads experimental mutual information input at channel capacity
exp_mi_input = np.loadtxt(expt_mi_input_path,delimiter=",")

# Loads single cell mutual information trajectories
mi_single_cells = np.loadtxt(mi_single_cell_path,delimiter=",")
# _____ Loading files END _____


# _____ Data Processing BEGIN _____
# Obtains the average and std dev of channel capacity input distribution across all instances of model
# obtained calculations
# Each instance being the result of 300 drawn cells
avg_cmi = np.average(cc_ceemi_input,axis = 0)
var_cmi = np.var(cc_ceemi_input,axis = 0)**.5
# _____ Data Processing END _____


# _____ Figure generation BEGIN _____
# Specifies figure colors
mdl_color = user_arguments["mdl_color"]
expt_color = user_arguments["expt_color"]
lim_color = user_arguments["lim_color"]

# Spread percentiles
sp = user_arguments["sp"]

# specifies times to plot
times_to_plot = [0, 6, 12, 24, 45, 60, 90]

# generates figure
fig, (ax) = plt.subplots(1, 1)

# Inset gives the cc at 90 minutes
axins = ax.inset_axes([.65,0.17, .3, .18])
# Labels for doses used
x_tick_labels = ["0","10","15","20","50"]
bar_width = .4
x_ticks = np.arange(len(x_tick_labels))
axins.bar(x_ticks-bar_width/2, exp_mi_input,bar_width, color = expt_color, label = "Experimental CSAR")
axins.bar(x_ticks+bar_width/2,avg_cmi,bar_width , yerr=var_cmi, color=mdl_color, label = "CeeMI")
axins.errorbar(x_ticks+ bar_width/2, avg_cmi,var_cmi ,fmt = "none",capsize=5,color = "k")
axins.set_xticks(x_ticks)
axins.set_xticklabels(x_tick_labels)
axins.set_ylim([0,0.55])
axins.set_xlim([-0.5,4.5])
axins.set_ylabel("Prob.")
axins.set_xlabel("Input (pM)")
axins.spines["top"].set_visible(False)
axins.spines["right"].set_visible(False)

# plots max input information
ax.hlines(np.log2(len(x_tick_labels)), -10,100, color=lim_color,linestyles="--")

# plots mutual information of experimental CSAR
ax.plot(times_to_plot, exp_mi, color=expt_color, label="CSAR")

# plots cell state conditioned mutual information
ax.errorbar(times_to_plot, np.mean(ceemi, axis=0), np.std(ceemi, axis=0), color=mdl_color, alpha=1,
             label="CeeMI")

# plots cell spread
ax.fill_between(times_to_plot, np.percentile(mi_single_cells, sp[0], axis=0),
                 np.percentile(mi_single_cells, sp[1], axis=0), color=mdl_color, alpha=0.3,
                 label="SC Spread", edgecolor="none")
ax.legend(loc = "upper center",ncol = 3,frameon=False, bbox_to_anchor=(.5,1.12))
ax.set_ylabel("Information (Bits)")
ax.set_xlabel("Time (Minutes)")
ax.set_ylim([0,2.35])
ax.set_xlim([0,90])
ax.text(25,2.2,"Max Information of Inputs",color = lim_color)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# generates and saves plot in the tight layout
plt.tight_layout()
plt.savefig(output_path)
plt.show()
# _____ Figure generation END _____
