import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path, PureWindowsPath

# Generates plot for experimental dose/response curve

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
# Specifies path to file which specifies user specified arguments for mapping cell parameters to cell responses
dose_response_dir = Path("Data/IGFR/Moments/experimental_moments_population/experimental_single_dose_all_dose_15_interval.csv")

# Specifies output path
output_path = Path("Figures/Supplementary/EGFR_Expt_Dose_Resp/dose_response.png")
# _____ File path declarations END _____


# _____ Loading files BEGIN _____
# loads dose response data
data = np.loadtxt(dose_response_dir,delimiter=',')
# _____ Loading files END _____


# _____ Figure generation BEGIN _____
# numerical dose list
dose = [0.0001,0.0078125,0.015625,0.03125,0.0625,0.125,0.25,0.5,1,100]

# doses to plot histogram x values
dose_to_hist = [0.0001,0.015625,0.5]

# doses to plot histogram dictionary keys
dose_hist_key = ['egf_0.0_nM','egf_0.015625_nM','egf_0.5_nM']

# specifies colors of the plot
mean_color = "black"
dev_color = "#FA4616"
dist_color = "#0021A5"

# list of stds
std = np.zeros(len(dose))

# list of mean responses
mean_response = np.zeros(len(dose))

# list of inputs
input_list = list(data.item().keys())

# populates arrays of standard deviation and mean responses
count = 0
for key in data.item().keys():
	std[count] = (np.std(data.item()[key]))
	mean_response[count] = np.mean(data.item()[key])
	count+=1
	
# generates figure
fig, ax = plt.subplots(1)

# plots standard deviation around mean response
ax.fill_between(dose,mean_response+std,mean_response-std,alpha = 0.5,label = "$\sigma$",color = dev_color)

# plots mean response
ax.plot(dose,mean_response,color = mean_color)

# Adds dots at measured values
ax.scatter(dose,mean_response,color = mean_color)

# Sets up axis for histograms which need different x axis
ax2=ax.twiny()

# obtains placement of base of histograms
max_r = 1E4
min_r = 0
ax2.set_xlim([min_r,max_r])
ax.set_xlim([0.00004,150])
display_dose_range_log_10 = np.log10(150) - np.log10(0.00004)
adjusted_doses = []
for d_t_h in dose_to_hist:
	adjusted_doses.append((max_r-min_r)*(np.log10(d_t_h)-np.log10(0.00004))/display_dose_range_log_10)

# plots histograms
for key_ind in range(len(dose_hist_key)):
	ax2.hist(data.item()[dose_hist_key[key_ind]],orientation = "horizontal",bottom = adjusted_doses[key_ind],alpha = .5,color = dist_color)
	
# sets x scale
ax.set_xscale('log')

# adds legend
ax.legend(frameon=False)

# adds labels
ax.set_xlabel("Extracellular EGF (nM)")
ax.set_ylabel("Surface EGFR $\left(10^{5}\\right)$")

# sets ticks
ax.set_yticks([0,1E5,3E5,5E5])
ax.set_yticklabels([0,1,3,5])
ax.set_xticks([0.0001,0.01,1,100])
ax.set_xticklabels([0,"$10^{-2}$","$10^{0}$","$10^{2}$"])

# removes x ticks from top of plot
ax2.set_xticks([])

# adds "Count" text to histograms
ax.text(0.0001,0.46E5,"Count")
ax.text(.015,0.46E5,"Count")
ax.text(.5,0.055E5,"Count")

# removes spines
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)

# generates plot in tight layout
plt.tight_layout()
plt.savefig(output_path)
plt.savefig(f"{output_path}.svg",format = 'svg')
plt.show()
