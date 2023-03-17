import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import scipy.stats as st
import os
from pathlib import Path

# Generates plot for EGFR response with both moment ode and gillespie algorithm

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
# Specifies path to file which contains results of gillespie simulations
gillespie_results_path = Path("Data/EGFR/Gillespie_Results/cell_0.npy")

# Specifies path to file which contains means from ode solution
mdl_means_path = Path("Data/EGFR/means_and_variances/noise_factor_1/mu_226.csv")

# Specifies path to file which contains variances from ode solution
mdl_vars_path = Path("Data/EGFR/means_and_variances/noise_factor_1/var_226.csv")

# Specifies path to file which contains headers of gillespie simulation
g_header_path = Path("Data/EGFR/Gillespie_Results/cell_header.txt")

# Specifies output path
output_path = Path("Figures/Supplementary/EGFR_Gillespie")

# Name of output image file
output_name = "EGFR_Results"
# _____ File path declarations END _____


# _____ Loading files BEGIN _____
# loads results of gillespie simulations
gillespie_results = np.load(gillespie_results_path)

# loads means from ode solution
mdl_means = np.loadtxt(mdl_means_path,delimiter=',')[0,:]

# loads variances from ode solution
mdl_vars = np.loadtxt(mdl_vars_path,delimiter=',')[0,:]

with open(g_header_path,'r') as header:
	g_header = header.readlines()
# _____ Loading files END _____

# _____ Figure generation BEGIN _____
# Specifies ode color
ode_color = '#FF5733'
# Specifies ssa color
ssa_color = '#8BD3E6'
# dimensions of figure
fig_dim = [7,10]
fig,axs = plt.subplots(10,1,sharex=True,figsize=(fig_dim[0], fig_dim[1]))
for i in range(np.size(mdl_means)):
	
	ssa_histogram = axs[i].hist(gillespie_results[:,i],density='pdf',bins=20,color = ssa_color)
	label = f"{round(eval(g_header[i])['dose'],3)} ng/mL EGF"
	axs[i].text(0,.99, label, ha='left', va='top', transform=axs[i].transAxes)
	axs[i].spines["top"].set_visible(False)
	axs[i].spines["right"].set_visible(False)
	axs[i].spines["left"].set_visible(False)
	axs[i].axes.yaxis.set_visible(False)
	

mdl_scale = mdl_vars/mdl_means
mdl_shape = mdl_means/mdl_scale
mins = np.min(gillespie_results,axis=0)
maxs = np.max(gillespie_results,axis=0)


for i in range(np.size(mdl_means)):
	# Note both distributions are shifted by the Gillespie mean
	shift = np.mean(gillespie_results[:,i])
	x = np.linspace(mins[i], maxs[i], 10000)
	pdf = st.gamma.pdf(x,mdl_shape[i],scale = mdl_scale[i])
	ode_plot = axs[i].plot(x,pdf,label="ODE",color = ode_color,linewidth=1)
	if i == 0:
		rectangle = Rectangle((0,0), 1, 1, color=ssa_color,label='Gillespie')
		axs[i].legend(handles = [rectangle,ode_plot[0]],frameon=False,ncol=2,loc='upper center', bbox_to_anchor=(0, 1, 1, 1))
plt.xlabel("sEGFR at Steady State")
#plt.xticks(ticks = [-1000,0,1000],labels = ["$\mu$-1000","$\mu$","$\mu$+1000"])
plt.tight_layout()
plt.savefig(output_path/output_name)
plt.savefig(output_path/(output_name+".svg"),format = 'svg')
plt.show()
# _____ Figure generation END _____
