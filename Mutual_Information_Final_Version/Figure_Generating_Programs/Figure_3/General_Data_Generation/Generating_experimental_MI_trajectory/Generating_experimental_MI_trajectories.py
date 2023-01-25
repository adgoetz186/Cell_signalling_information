import numpy as np
import os
import shutil
import Mutual_Information_Final_Version.functions.data_processing_functions.Conditioning_Functions as cf
import Mutual_Information_Final_Version.functions.data_processing_functions.Process_Raw_Cell_Conditioning as prc
import Mutual_Information_Final_Version.functions.mutual_information_functions.MI_Calculation_Functions as mi
import Mutual_Information_Final_Version.functions.file_processing_functions.user_input_functions as uif
import Mutual_Information_Final_Version.functions.file_processing_functions.make_folder as mk

# This program generates the experimental MI of the CSAR trajectory for the igfr/akt/foxo pathway

# _____ Load arguments BEGIN _____
# Specifies path to file which specifies user specified arguments
user_input_file_path = "/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Figure_Generating_Programs/Figure_3/General_Data_Generation/Generating_experimental_MI_trajectory/args.txt"
user_arguments = uif.read_file_to_arg_dict(user_input_file_path)
# _____ Load arguments END _____

# _____ File path declarations BEGIN _____
# Specifies path to experimental moment file
moments_file_path = user_arguments["moments_file_path"]

# Specifies path to headers for experimental moment file
cell_moments_header_path = user_arguments["cell_moments_header_path"]

# Specifies path to output directory
output_dir = user_arguments["output_dir"]

# Specifies name of input distribution file
input_dist_filename = user_arguments["input_dist_filename"]

# Specifies name of information performance file
info_perf_filename = user_arguments["info_perf_filename"]
# _____ File path declarations END _____


# _____ Loading files BEGIN _____
# Loads experimental moments
moments = np.loadtxt(moments_file_path,delimiter=",")
# _____ Loading files END _____


# _____ Main code BEGIN _____
# specifies single cell response type and discretization parameter
discretization_parameter = user_arguments["discretization_parameter"]
distribution_type = user_arguments["distribution_type"]

# specifies number of doses
dose_count = user_arguments["dose_count"]

# provides list of times to plot
time_list = user_arguments["time_list"]

# provides time to use to save the cc input distribution for plotting
time_for_cc_input = user_arguments["time_for_cc_input"]

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
