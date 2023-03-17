import numpy as np
import os
import shutil
from pathlib import Path, PureWindowsPath
import Mutual_Information_Main.functions.data_processing_functions.Conditioning_Functions as cf
import Mutual_Information_Main.functions.data_processing_functions.Process_Raw_Cell_Conditioning as prc
import Mutual_Information_Main.functions.mutual_information_functions.MI_Calculation_Functions as mi
import Mutual_Information_Main.functions.file_processing_functions.make_folder as mk

# This program generates the experimental I trajectory for the igfr/akt/foxo pathway

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
# Specifies path to experimental moment file
moments_file_path = Path("Data/IGFR/Moments/experimental_moments_population/experimental_single_dose.csv")

# Specifies path to headers for experimental moment file
cell_moments_header_path = Path("Data/IGFR/Moments/experimental_moments_population/experimental_single_dose_header.txt")

# Specifies path to output directory
output_dir = Path("Data/IGFR/Information_Values/Experimental_MI_CSAR/")

# Specifies name of input distribution file
input_dist_filename = "exp_cc_mi"

# Specifies name of information performance file
info_perf_filename = "exp_mi"
# _____ File path declarations END _____


# _____ Loading files BEGIN _____
# Loads experimental moments
moments = np.loadtxt(moments_file_path,delimiter=",")
# _____ Loading files END _____


# _____ Main code BEGIN _____
# specifies number of doses
dose_count = 5

# provides list of times to plot
time_list = [0, 6, 12, 24, 45, 60, 90]

# provides time to use to save the cc input distribution for plotting
time_for_cc_input = 90

# initialize the array to store information values
mi_all_cell_array = np.zeros(len(time_list))
for time_index in range(len(time_list)):
	# Intermediate folder locations/names
	exp_MI_raw_data_folder_name = f"{user_arguments['exp_MI_raw_data_folder_prefix']}_time_{time_list[time_index]}"
	exp_MI_raw_data_folder_location = user_arguments["exp_MI_raw_data_folder_location"]
	exp_MI_crm_data_folder_location = user_arguments["exp_MI_crm_data_folder_location"]
	exp_MI_crm_data_folder_name = f"{user_arguments['exp_MI_raw_data_folder_prefix']}_time_{time_list[time_index]}"
	
	# if path does not exists, makes new folders
	if not os.path.exists(exp_MI_raw_data_folder_location):
		mk.make_folder(f"{exp_MI_raw_data_folder_location}")
	if not os.path.exists(exp_MI_crm_data_folder_location):
		mk.make_folder(f"{exp_MI_crm_data_folder_location}")
	
	# creates raw cell conditioning folder and fills it with raw cell files (moments)
	cf.create_dir_of_cell_data_files({}, moments, f"{exp_MI_raw_data_folder_location}/{exp_MI_raw_data_folder_name}", cell_data_column_header_name_and_location=cell_moments_header_path)
	
	# creates and populates folder of conditional response matrices
	mk.make_folder(f"{exp_MI_crm_data_folder_location}/{exp_MI_crm_data_folder_name}")
	for cell_file_name in os.listdir(f"{exp_MI_raw_data_folder_location}/{exp_MI_raw_data_folder_name}"):
		prc.moment_file_to_crm_file(f"{exp_MI_raw_data_folder_location}/{exp_MI_raw_data_folder_name}/{cell_file_name}", f"{exp_MI_crm_data_folder_location}/{exp_MI_crm_data_folder_name}/{cell_file_name}",lock_parameter={"time": f"{time_list[time_index]}_min"},assumed_distribution=distribution_type, discretization_parameter=discretization_parameter)
	
	# obtains channel capacity for Cee-MI
	mi_all_cell, mi_cc_in = mi.pcmi_at_cc(f"{exp_MI_crm_data_folder_location}/{exp_MI_crm_data_folder_name}", return_all_cell_values=True)
	mi_all_cell_array[time_index] = mi_all_cell
	
	# stores predefined cc input array
	if time_list[time_index] == time_for_cc_input:
		input_dist = mi_cc_in

# creates and fills output folder
if not os.path.exists(output_dir):
	mk.make_folder(f"{output_dir}")
np.savetxt(f"{output_dir}/{input_dist_filename}", input_dist, delimiter=",")
np.savetxt(f"{output_dir}/{info_perf_filename}", mi_all_cell_array , delimiter=",")
# _____ Main code End _____
