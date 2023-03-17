import numpy as np
import json
import scipy.stats as st
import matplotlib.pyplot as plt
import os
from pathlib import Path

# Generates plot for IGFR response with both moment ode and gillespie algorithm

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
gillespie_results_path = Path("Data/IGFR/Gillespie_Results/cell_0.npy")

# Specifies path to file which contains moments from ode solution
mdl_moments_path = Path("Data/IGFR/Moments/Model_Moments/3_dose_response/moments_162_3_dose_full_time.npy")

# Specifies path to file which contains headers of gillespie simulation
g_header_path = Path("Data/IGFR/Gillespie_Results/cell_header.txt")

# Specifies output path
output_path = Path("Figures/Supplementary/IGFR_Gillespie")

# Specifies output name
output_name = "IGFR_Results"
# _____ File path declarations END _____


# _____ Loading files BEGIN _____
# loads results of gillespie simulations
gillespie_results = np.load(gillespie_results_path)

# loads means from ode solution
moment_results = np.load(mdl_moments_path)[0]

with open(g_header_path,'r') as header:
	g_header = header.readlines()
# _____ Loading files END _____

# _____ Figure generation BEGIN _____
g_header = [i.replace("'",'"') for i in g_header]
g_header = [json.loads(i) for i in g_header]
print(g_header)
print(gillespie_results)
print(g_header[:31])
dose_1 = gillespie_results[:,:31]
dose_2 = gillespie_results[:,31:62]
dose_3 = gillespie_results[:,62:]
print(dose_1)
time_list = list(set([i['time'] for i in g_header]))
time_list.sort()
print(time_list)



times = np.arange(0,90+3,3)*60
print(times)
print(np.shape(moment_results))
print(moment_results[93:124])
print(moment_results[:31])

m_dose_1_mean = moment_results[:31]
m_dose_1_var = moment_results[93:124]-m_dose_1_mean**2
m_dose_1_scale = m_dose_1_var/m_dose_1_mean
m_dose_1_shape = m_dose_1_mean/m_dose_1_scale



m_dose_2_mean = moment_results[31:62]
m_dose_2_var = moment_results[124:155]-m_dose_2_mean**2
m_dose_2_scale = m_dose_2_var/m_dose_2_mean
m_dose_2_shape = m_dose_2_mean/m_dose_2_scale

m_dose_3_mean = moment_results[62:93]
m_dose_3_var = moment_results[155:186]-m_dose_3_mean**2
m_dose_3_scale = m_dose_3_var/m_dose_3_mean
m_dose_3_shape = m_dose_3_mean/m_dose_3_scale

# dimensions of figure
fig_dim = [7,10]

fig, axs = plt.subplots(7,1,sharex=True,figsize=(fig_dim[0], fig_dim[1]))

time_indicies = [0,2,4,8,15,20,30]

for i in range(len(time_indicies)):
	print(i)
	print(m_dose_1_shape)
	print(m_dose_1_shape[time_indicies[i]])
	x = np.linspace(np.min(dose_1[:,time_indicies[i]]),np.max(dose_1[:,time_indicies[i]]),10000)
	gamma_dist = st.gamma.pdf(x,m_dose_1_shape[time_indicies[i]],scale = m_dose_1_scale[time_indicies[i]])
	axs[i].plot(x,gamma_dist,color = (0,.7,1),zorder=10,linewidth = 3)
	
	x = np.linspace(np.min(dose_2[:, time_indicies[i]]), np.max(dose_2[:, time_indicies[i]]), 10000)
	gamma_dist = st.gamma.pdf(x, m_dose_2_shape[time_indicies[i]], scale=m_dose_2_scale[time_indicies[i]])
	axs[i].plot(x, gamma_dist,color = (1,.5,.5),zorder=11,linewidth = 3)
	
	
	x = np.linspace(np.min(dose_3[:, time_indicies[i]]), np.max(dose_3[:, time_indicies[i]]), 10000)
	gamma_dist = st.gamma.pdf(x, m_dose_3_shape[time_indicies[i]], scale=m_dose_3_scale[time_indicies[i]])
	axs[i].plot(x, gamma_dist,color = (.1,1,0),zorder=12,linewidth = 3)
	
	axs[i].hist(dose_1[:, time_indicies[i]], density='pdf',
	            bins=np.arange(np.min(dose_1[:, time_indicies[i]]), np.max(dose_1[:, time_indicies[i]]) + 1),
	            color="blue",zorder=10,label="0 pM IGF")
	axs[i].hist(dose_2[:, time_indicies[i]], density='pdf',
	            bins=np.arange(np.min(dose_2[:, time_indicies[i]]), np.max(dose_2[:, time_indicies[i]]) + 1),
	            color="orange",zorder=11,label="20 pM IGF")
	axs[i].hist(dose_3[:, time_indicies[i]], density='pdf',
	            bins=np.arange(np.min(dose_3[:, time_indicies[i]]), np.max(dose_3[:, time_indicies[i]]) + 1),
	            color="green",zorder=12,label="50 pM IGF")
	axs[i].spines["top"].set_visible(False)
	axs[i].spines["right"].set_visible(False)
	axs[i].spines["left"].set_visible(False)
	axs[i].axes.yaxis.set_visible(False)
	axs[i].text(0, .99, f"{time_indicies[i]*3} Min", ha='left', va='top', transform=axs[i].transAxes)
	if i == 0:
		#rectangle = Rectangle((0, 0), 1, 1, color=ssa_color, label='SSA')
		axs[i].legend(frameon=False, ncol=3, loc='upper center',
		              bbox_to_anchor=(0, 1, 1, 1))
plt.xlabel('Nuclear Foxo Levels')
plt.tight_layout()
plt.savefig(output_path / output_name)
plt.savefig(output_path / (output_name+".svg"),format = 'svg')
plt.show()
#ax1.plot(time_list,np.mean(dose_1,axis=0))
#ax1.fill_between(time_list,np.mean(dose_1,axis=0)+np.std(dose_1,axis=0),np.mean(dose_1,axis=0)-np.std(dose_1,axis=0),alpha=0.4)
#ax1.plot(time_list,np.mean(dose_2,axis=0))
#ax1.fill_between(time_list,np.mean(dose_2,axis=0)+np.std(dose_2,axis=0),np.mean(dose_2,axis=0)-np.std(dose_2,axis=0),alpha=0.4)
#ax1.plot(time_list,np.mean(dose_3,axis=0))
#ax1.fill_between(time_list,np.mean(dose_3,axis=0)+np.std(dose_3,axis=0),np.mean(dose_3,axis=0)-np.std(dose_3,axis=0),alpha=0.4)
#ax1.spines["top"].set_visible(False)
#ax1.spines["right"].set_visible(False)


#ax2.plot(times,m_dose_1_mean)
#ax2.fill_between(times,m_dose_1_mean+m_dose_1_var**.5,m_dose_1_mean-m_dose_1_var**.5,alpha=0.4)
#ax2.plot(times,m_dose_2_mean)
#ax2.fill_between(times,m_dose_2_mean+m_dose_2_var**.5,m_dose_2_mean-m_dose_2_var**.5,alpha=0.4)
#ax2.plot(times,m_dose_3_mean)
#ax2.fill_between(times,m_dose_3_mean+m_dose_3_var**.5,m_dose_3_mean-m_dose_3_var**.5,alpha=0.4)
#ax2.spines["top"].set_visible(False)
#ax2.spines["right"].set_visible(False)

#plt.show()
# _____ Figure generation END _____