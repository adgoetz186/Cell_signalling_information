import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import os
from pathlib import Path, PureWindowsPath

# Generates plot of EGFR signaling performance with difference levels of noise and MI of cell state agnostic response

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
# location of Gillespie simulation results
ssa_path = Path("Data/Toy_Model/Gillespie_Results/cell_0.npy")
# Specifies path to file which contains headers of gillespie simulation
g_header_path = Path("Data/Toy_Model/Gillespie_Results/cell_header.txt")
# location of output image
output_path = Path("Figures/Supplementary/Toy_Gillespie")
# name of output image
output_name = "Toy_Gillespie"
# _____ File path declarations END _____

# _____ Loading files BEGIN _____
# loads the results of the ssa
ssa = np.load(ssa_path)

with open(g_header_path,'r') as header:
	g_header = header.readlines()
# _____ Loading files End _____

# _____ Figure generation BEGIN _____
# Values should be the same as the ones given in the gillespie simulation
# average degredation rate, k_deg
deg_mean = 5
# average initial receptor count value, r0
r0_mean = 500
# rate of binding, k_bind
k_bind = 1
# rate of unbinding, k_unbind
k_unbind = 10
# Specifies ode color
ode_color = '#FF5733'
# Specifies ssa color
ssa_color = '#8BD3E6'
# dimensions of figure
fig_dim = [4.2,8]
fig, axs = plt.subplots(9, 1, sharex=True, figsize=(fig_dim[0], fig_dim[1]))
for subplot_number in range(9):
	i = subplot_number * 3
	ssa_histogram = axs[subplot_number].hist(ssa[:, i], density='pdf', bins=20, color=ssa_color)
	label = f"{round(eval(g_header[i])['dose'], 1)} a.u. L"
	if subplot_number < 5:
		axs[subplot_number].text(1, 1, label, ha='right', va='top', transform=axs[subplot_number].transAxes)
	else:
		axs[subplot_number].text(0, 1, label, ha='left', va='top', transform=axs[subplot_number].transAxes)
	axs[subplot_number].spines["top"].set_visible(False)
	axs[subplot_number].spines["right"].set_visible(False)
	axs[subplot_number].spines["left"].set_visible(False)
	axs[subplot_number].axes.yaxis.set_visible(False)

ulist = np.array([eval(i)['dose'] for i in g_header])
mdl_rates = ulist*k_bind*r0_mean/(ulist*k_bind+deg_mean+k_unbind)
mins = np.min(ssa, axis=0)
maxs = np.max(ssa, axis=0)

for subplot_number in range(9):
	i = subplot_number * 3
	x = np.arange(mins[i], maxs[i])
	pmfv = st.poisson.pmf(x, mdl_rates[i])
	ode_plot = axs[subplot_number].plot(x, pmfv, label="Poisson", color=ode_color, linewidth=2)
	if i == 0:
		rectangle = Rectangle((0, 0), 1, 1, color=ssa_color, label='Gillespie')
		axs[subplot_number].legend(handles=[rectangle, ode_plot[0]], frameon=False, ncol=1, loc='upper center',
		              bbox_to_anchor=(0, 1, 1, 1))
plt.xlabel("Bound Receptor Count (Steady State)")
#plt.xticks(ticks=[-100, 0, 100], labels=["$\mu$-100", "$\mu$", "$\mu$+100"])
plt.tight_layout()
plt.savefig(output_path / output_name)
plt.savefig(output_path / (output_name + ".svg"), format='svg')
plt.show()
# _____ Figure generation END _____
