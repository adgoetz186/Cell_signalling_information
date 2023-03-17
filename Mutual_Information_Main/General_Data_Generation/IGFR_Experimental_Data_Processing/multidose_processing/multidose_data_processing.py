import os
import numpy as np
import copy as cp
from pathlib import Path, PureWindowsPath

# Processes the raw multi-dose florescence data provided by Heiser Lab
# https://pubmed.ncbi.nlm.nih.gov/31838146/

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
# location for multidose file as provided by Heiser Lab
raw_mutlidose_file = Path("Data/IGFR/experimental_data/Raw/MultiDose.csv")
# location to store processed multidose file
output_mutlidose_file = "Data/IGFR/experimental_data/adjusted_multidose/multi_dose"
# _____ File path declarations END _____

# _____ Loading files BEGIN _____
mtd_data = np.loadtxt(raw_mutlidose_file,delimiter=",")
# _____ Loading files End _____

# _____ Main code BEGIN _____
# Background measurements for each cell to remove
background = mtd_data[:,0]

# Time Values
times = mtd_data[:,1]

# Measured flourescence trajectories
trajectories = mtd_data[:,2:]

initial_n_foxo = 2/3*710

# Removes measured background values
background_removed_trajectories = trajectories - np.reshape(background,(-1,1))

# converts from a.u. to average protein count
background_removed_trajectories *= initial_n_foxo/np.average(background_removed_trajectories[0])

# rearranges and saves array
background_removed_trajectories = np.transpose(background_removed_trajectories)
background_removed_trajectories = np.vstack((times,background_removed_trajectories))
np.save(output_mutlidose_file,background_removed_trajectories)
# _____ Main code END _____
