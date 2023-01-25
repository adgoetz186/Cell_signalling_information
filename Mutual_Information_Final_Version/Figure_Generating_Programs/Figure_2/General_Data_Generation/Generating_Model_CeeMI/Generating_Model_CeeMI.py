import numpy as np
import os
import Mutual_Information_Final_Version.functions.data_processing_functions.Conditioning_Functions as cf
import Mutual_Information_Final_Version.functions.data_processing_functions.Process_Raw_Cell_Conditioning as prc
import Mutual_Information_Final_Version.functions.mutual_information_functions.MI_Calculation_Functions as mi
import Mutual_Information_Final_Version.functions.file_processing_functions.user_input_functions as uif
import Mutual_Information_Final_Version.functions.file_processing_functions.make_folder as mk

# Generates cell state dependent mutual information values for EGFR model fit data

# _____ Load arguments BEGIN _____
# Specifies path to file which specifies user specified arguments
user_input_file_path = "/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Figure_Generating_Programs/Figure_2/General_Data_Generation/Generating_Model_CeeMI/noise_1_args.txt"
user_arguments = uif.read_file_to_arg_dict(user_input_file_path)
# _____ Load arguments END _____

# _____ File path declarations BEGIN _____
# Specifies path to experimental moment file
means_dir = user_arguments["mean_dir"]
var_dir = user_arguments["var_dir"]

# Output file location
output_information_dir = user_arguments["output_information_dir"]

# Intermediate folder locations/names
CeeMI_raw_data_folder_location = user_arguments["CeeMI_raw_data_folder_location"]
CeeMI_crm_data_folder_location = user_arguments["CeeMI_crm_data_folder_location"]
CeeMI_raw_data_folder = user_arguments["CeeMI_raw_data_folder"]
CeeMI_crm_data_folder = user_arguments["CeeMI_crm_data_folder"]
# _____ File path declarations END _____

# _____ Loading files BEGIN _____
mean_array = np.loadtxt(means_dir,delimiter=",")
var_array = np.loadtxt(var_dir,delimiter=",")
# _____ Loading files End _____

# _____ Main code BEGIN _____
# specifies the number of cells to simulate
cell_count = user_arguments["cell_count"]

# combines moment arrays
moment_array = np.hstack((mean_array,mean_array**2+var_array))[:cell_count]

# creates raw cell conditioning folder and fills it with raw cell files (moments)
cf.create_dir_of_cell_data_files({}, moment_array,f"{CeeMI_raw_data_folder_location}/{CeeMI_raw_data_folder}",all_cells_are_conditioned_on=True)

# creates and populates folder of conditional response matrices
mk.make_folder(f"{CeeMI_crm_data_folder_location}/{CeeMI_crm_data_folder}")
for cell_file_name in os.listdir(f"{CeeMI_raw_data_folder_location}/{CeeMI_raw_data_folder}"):
	prc.moment_file_to_crm_file(f"{CeeMI_raw_data_folder_location}/{CeeMI_raw_data_folder}/{cell_file_name}",f"{CeeMI_crm_data_folder_location}/{CeeMI_crm_data_folder}/{cell_file_name}",assumed_distribution=user_arguments["distribution_type"],discretization_parameter=user_arguments["discretization_parameter"])

# obtains channel capacity for Cee-MI
cmi_all_cell, cmi_cc_in = mi.pcmi_at_cc(f"{CeeMI_crm_data_folder_location}/{CeeMI_crm_data_folder}", return_all_cell_values=True)

# saves output
if not os.path.exists(output_information_dir):
	os.mkdir(output_information_dir)
np.savetxt(f"{output_information_dir}/single_cell_MI.csv",cmi_all_cell,delimiter=",")
np.savetxt(f"{output_information_dir}/cc_input_distribution",cmi_cc_in,delimiter=",")
# _____ Main code END _____
