import random
import numpy as np
import os
import shutil
from pathlib import Path, PureWindowsPath
import Mutual_Information_Main.functions.data_processing_functions.Process_Raw_Cell_Conditioning as prc
import Mutual_Information_Main.functions.mutual_information_functions.MI_Calculation_Functions as mi

# This program generates the model I(\theta) values for the igfr/akt/foxo pathway

# _____ Setting the CWD to be Mutual_Information_Main BEGIN _____
# Cell_signaling_information path here
path_to_CSI = ""
if path_to_CSI == "":
	try:
		# Obtains the location of the Cell_signaling_information folder if it is in the cwd parents
		path_to_CSI = Path.cwd().parents[
			[Path.cwd().parents[i].parts[-1] for i in range(len(Path.cwd().parts) - 1)].index(
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
path_to_MIM = path_to_CSI / "Mutual_Information_Main"
os.chdir(path_to_MIM)
# _____ Setting the CWD to be Mutual_Information_Main END _____

# _____ File path declarations BEGIN _____
# Specifies path to drawn model cell parameters
params_file_path = Path("Data/IGFR/Model_Params/params_162/params_162.npy")
# Specifies path to drawn model cell moments
moments_file_path = Path("Data/IGFR/Moments/Model_Moments/4_dose_response/moments_162_4_dose_60_90_min.npy")
# Specifies path to headers for model parameter file
cell_params_header_path = Path("Data/IGFR/Model_Params/params_162/params_162_header_structure.txt")
# Specifies path to headers for model moment file
cell_moments_header_path = Path(
	"Data/IGFR/Moments/Model_Moments/4_dose_response/moments_162_4_dose_60_90_min_header_structure.csv")
# Specifies path to output directory
output_dir = Path("Data/IGFR/Information_Values/Model_CeeMI_Step_Dose")
# _____ File path declarations END _____

# _____ Loading files BEGIN _____
# Loads parameters of drawn cells
params = np.load(params_file_path)

# Loads moments of drawn cells
moments = np.load(moments_file_path)
# _____ Loading files END _____


# _____ Main code BEGIN _____
even_tpl = tuple([2 * i for i in range(np.shape(moments)[1] // 2)])
odd_tpl = tuple([2 * i + 1 for i in range(np.shape(moments)[1] // 2)])

# Averages 60 and 90 minute moments
moments = (moments[:, even_tpl] + moments[:, odd_tpl]) / 2

all_parameter_list = list(np.loadtxt(cell_params_header_path, dtype=str))
cell_index_list = list(range(np.shape(moments)[0]))

# number of cells (50000)
number_of_cells = 50000
number_of_inputs = 4
random.shuffle(cell_index_list)

selected_cells = cell_index_list[:number_of_cells]
params_select = params[selected_cells, :]
single_cell_mi_array = np.zeros_like(moments)
single_cell_mi_params = np.zeros_like(params)

cc_ceemi_array = np.zeros((1, number_of_inputs))

moments_select = moments[selected_cells, :]

list_of_param_values = [10 ** params_select[:, i] for i in range(np.shape(params_select)[1])]
param_dict = dict(zip(all_parameter_list, list_of_param_values))

moments_list = prc.create_list_of_cell_data_files(param_dict, moments_select, all_cells_are_conditioned_on=True)

crm_list = prc.moment_list_to_crm_list(moments_list)

cmi_all_cell, cc_in = mi.pcmi_at_cc_from_crm_list(crm_list, return_all_cell_values=True)

if not os.path.exists(output_dir):
	os.mkdir(output_dir)
np.save(output_dir / "mdl_moments_used", moments_select)
np.save(output_dir / "mdl_params_used", params_select)
np.save(output_dir / "cc_in", cc_in)
np.save(output_dir / "single_cell_mi_array", cmi_all_cell)
# _____ Main code END _____
