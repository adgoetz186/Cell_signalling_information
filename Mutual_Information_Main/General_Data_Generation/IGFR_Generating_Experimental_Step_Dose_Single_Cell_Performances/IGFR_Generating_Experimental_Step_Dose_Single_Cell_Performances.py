import matplotlib.pyplot as plt
import numpy as np
import os
import shutil
import scipy.stats as st
from pathlib import Path, PureWindowsPath
import Mutual_Information_Main.functions.data_processing_functions.Process_Raw_Cell_Conditioning as prc
import Mutual_Information_Main.functions.mutual_information_functions.MI_Calculation_Functions as mi

# This program generates the experimental I(\theta) values for the igfr/akt/foxo pathway in figure 2

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
# Specifies path to drawn model cell moments
moments_file_path = Path("Data/IGFR/Moments/experimental_moments_single_cell/multidose_moments.npy")
# Specifies path to headers for model moment file
cell_moments_header_path = Path("Data/IGFR/Moments/experimental_moments_single_cell/multidose_moments_header.txt")
# Specifies path to output directory
output_dir = Path("Data/IGFR/Information_Values/Experimental_Step_Dose_Single_Cell_Performances")
# Specifies output file name
output_file_name = "Expt_CeeMI_Step_Dose.csv"
# _____ File path declarations END _____
moments = np.load(moments_file_path)

# parameter_dict_fi = {"foxo": foxo, "igfr": igfr}

moment_list = prc.create_list_of_cell_data_files({}, moments, all_cells_are_conditioned_on=True)

crm_list = prc.moment_list_to_crm_list(moment_list)

# obtains optimal ICee
cmi_all_cell, cmi_cc_in = mi.pcmi_at_cc_from_crm_list(crm_list, return_all_cell_values=True)

# saves output
if not os.path.exists(output_dir):
	os.mkdir(output_dir)
np.save(output_dir/"single_cell_MI",cmi_all_cell)
np.save(output_dir/"cc_input_distribution",cmi_cc_in)
