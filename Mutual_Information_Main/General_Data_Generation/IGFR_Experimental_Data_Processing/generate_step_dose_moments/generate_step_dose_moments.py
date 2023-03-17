import os
import numpy as np
import copy as cp
from pathlib import Path, PureWindowsPath

# Obtains first 2 moments of step dose data at 4 steady states

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
# location for processed step dose data
processed_multidose_file = Path("Data/IGFR/experimental_data/adjusted_multidose/multi_dose.npy")
# location to store experimental moments
output_moment_file = Path("Data/IGFR/Moments/experimental_moments_single_cell/multidose_moments")
# _____ File path declarations END _____

# _____ Loading files BEGIN _____
mtd_data = np.load(processed_multidose_file)
# _____ Loading files End _____

# _____ Main code BEGIN _____
moment_file = np.zeros((np.shape(mtd_data)[0]-1,8))
moment_file[:,0] = np.average(mtd_data[1:,0:11],axis=1)
moment_file[:,1] = np.average(mtd_data[1:,30:41],axis=1)
moment_file[:,2] = np.average(mtd_data[1:,60:71],axis=1)
moment_file[:,3] = np.average(mtd_data[1:,90:101],axis=1)
moment_file[:,4] = np.var(mtd_data[1:,0:11],axis=1) + np.average(mtd_data[1:,0:11],axis=1)**2
moment_file[:,5] = np.var(mtd_data[1:,30:41],axis=1) + np.average(mtd_data[1:,30:41],axis=1)**2
moment_file[:,6] = np.var(mtd_data[1:,60:71],axis=1) + np.average(mtd_data[1:,60:71],axis=1)**2
moment_file[:,7] = np.var(mtd_data[1:,90:101],axis=1) + np.average(mtd_data[1:,90:101],axis=1)**2
np.save(output_moment_file,moment_file)
# _____ Main code END _____
