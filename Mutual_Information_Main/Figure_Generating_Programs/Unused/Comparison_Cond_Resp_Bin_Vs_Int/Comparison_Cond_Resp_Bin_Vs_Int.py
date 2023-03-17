import numpy as np
import time
import os
from pathlib import Path, PureWindowsPath
import scipy.stats as st
import matplotlib.pyplot as plt
from Mutual_Information_Main.functions.toy_model_functions import Toy_Model_Functions as tmf
import pandas as pd

# Generates schematic showing cell state agnostic responses and single cell responses

plt.rcParams.update({'font.size': 15})
plt.rcParams["hatch.linewidth"] = 4

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
# location of output folder for supplementary figures
output_dir = "Figures/Supplementary/Comparison_Cond_Resp_Bin_Vs_Int/"

# location of single cell egfr binned MI values
sc_egfr_bin_path = "Data/EGFR/moment_comparison_integration_methods/sc_binning.npy"

# location of single cell egfr integration MI values
sc_egfr_int_path = "Data/EGFR/moment_comparison_integration_methods/sc_integration.npy"

# location of single cell igfr binned MI values
sc_igfr_bin_path = "Data/IGFR/moment_comparison_integration_methods/sc_binning.npy"

# location of single cell igfr integration MI values
sc_igfr_int_path = "Data/IGFR/moment_comparison_integration_methods/sc_integration.npy"

# location of population egfr binned MI values
pop_egfr_bin_path = "Data/EGFR/moment_comparison_integration_methods/pop_binning.npy"

# location of population egfr integration MI values
pop_egfr_int_path = "Data/EGFR/moment_comparison_integration_methods/pop_integration.npy"

# location of population igfr binned MI values
pop_igfr_bin_path = "Data/IGFR/moment_comparison_integration_methods/pop_binning.npy"

# location of population igfr integration MI values
pop_igfr_int_path = "Data/IGFR/moment_comparison_integration_methods/pop_integration.npy"
# _____ File path declarations END _____

# _____ Loading files BEGIN _____
# Loads results of single cell egfr binning method
sc_egfr_bin = np.load(sc_egfr_bin_path)

# Loads results of single cell egfr integration method
sc_egfr_int = np.load(sc_egfr_int_path)

# Loads results of single cell igfr binning method
sc_igfr_bin = np.load(sc_igfr_bin_path)

# Loads results of single cell igfr integration method
sc_igfr_int = np.load(sc_igfr_int_path)

# Loads results of population egfr binning method
pop_egfr_bin = np.load(pop_egfr_bin_path)

# Loads results of population egfr integration method
pop_egfr_int = np.load(pop_egfr_int_path)

# Loads results of population igfr binning method
pop_igfr_bin = np.load(pop_igfr_bin_path)

# Loads results of population igfr integration method
pop_igfr_int = np.load(pop_igfr_int_path)
# _____ Loading files END _____

# _____ Figure generation BEGIN _____
fig,axs = plt.subplots(1,2,sharey=True,figsize=[8,4])
print(np.log10(np.abs(pop_egfr_bin-pop_egfr_int)/pop_egfr_int))
axs[0].hist(np.abs(sc_egfr_bin-sc_egfr_int)/sc_egfr_int*100,weights = np.ones_like(sc_egfr_bin)/np.size(sc_egfr_bin),bins=25,log=True)
print(np.abs(pop_igfr_bin-pop_igfr_int)/pop_igfr_int)
axs[1].hist(np.abs(sc_igfr_bin-sc_igfr_int)/sc_igfr_int*100,weights = np.ones_like(sc_igfr_bin)/np.size(sc_igfr_bin),bins=25,log=True)
axs[0].set_xticks([0,0.005,0.01])
axs[0].set_ylabel("Probability")
axs[0].spines["top"].set_visible(False)
axs[0].spines["right"].set_visible(False)
axs[1].spines["top"].set_visible(False)
axs[1].spines["right"].set_visible(False)
axs[0].set_xlabel("% Error EGFR")
axs[1].set_xlabel("% Error IGFR-Akt-FoxO")
plt.tight_layout()
plt.savefig(output_dir+"error_between_methods")
plt.savefig(output_dir+"error_between_methods.svg",format='svg')
plt.show()

# _____ Figure generation BEGIN _____
