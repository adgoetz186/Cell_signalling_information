import numpy as np
import os
from pathlib import Path, PureWindowsPath
import Mutual_Information_Main.functions.data_processing_functions.Process_Raw_Cell_Conditioning as prc
import Mutual_Information_Main.functions.mutual_information_functions.MI_Calculation_Functions as mi

# Generates mutual information value of cell state agnostic response from experimental data for EGFR

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
# Specifies path to experimental moment files
expt_mom_1 = Path("Data/EGFR/Experimental_Moments/EGFR_mean_10Conditions_20201116.npy")
expt_mom_2 = Path("Data/EGFR/Experimental_Moments/EGFR_2ndMomentMean_10Conditions_20201116.npy")

# Output file location
output_information_dir = Path("Data/EGFR/Information_Values/expt_EGFR")
# _____ File path declarations END _____


# _____ Loading files BEGIN _____
# loads moment files
expt_mom_1_array = np.load(expt_mom_1)
expt_mom_2_array = np.load(expt_mom_2)
# _____ Loading files End _____

# _____ Main code BEGIN _____
# conversion factor from a.u. to receptor count
au_cf = 0.00122

# converts moments from a.u. into appropriate units
expt_mom_1_array/=au_cf
expt_mom_2_array/=au_cf**2

# combines moment arrays
moment_array = np.hstack((expt_mom_1_array,expt_mom_2_array))

# creates list of moments which define each channel (cell state averaged here)
moment_list = prc.create_list_of_cell_data_files({}, moment_array)

# creates list of conditional response matrices for each channel
crm_list = prc.moment_list_to_crm_list(moment_list)

# obtains channel capacity for MI of CSAR
cmi_all_cell, cmi_cc_in = mi.pcmi_at_cc_from_crm_list(crm_list, return_all_cell_values=True)

# saves output
if not os.path.exists(output_information_dir):
	os.mkdir(output_information_dir)
np.save(output_information_dir/"single_cell_MI",cmi_all_cell)
np.save(output_information_dir/"cc_input_distribution",cmi_cc_in)
# _____ Main code End _____
