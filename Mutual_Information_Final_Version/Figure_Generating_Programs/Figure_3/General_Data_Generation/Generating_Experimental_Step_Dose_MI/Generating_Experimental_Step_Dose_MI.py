import matplotlib.pyplot as plt
import numpy as np
import os
import shutil
import scipy.stats as st
import Mutual_Information_Final_Version.functions.data_processing_functions.Conditioning_Functions as cf
import Mutual_Information_Final_Version.functions.data_processing_functions.Process_Raw_Cell_Conditioning as prc
import Mutual_Information_Final_Version.functions.file_processing_functions.user_input_functions as uif
import Mutual_Information_Final_Version.functions.mutual_information_functions.MI_Calculation_Functions as mi

# This program generates the model MI and CMI trajectories for the igfr/akt/foxo pathway in figure 2

# _____ Load arguments BEGIN _____
# Specifies path to file which specifies user specified arguments
user_input_file_path = "/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Figure_Generating_Programs/Figure_3/General_Data_Generation/Generating_Experimental_Step_Dose_MI/args.txt"
user_arguments = uif.read_file_to_arg_dict(user_input_file_path)
# _____ Load arguments END _____


# _____ File path declarations BEGIN _____
# Specifies path to drawn model cell moments
moments_file_path = user_arguments["moments_file_path"]
# Specifies path to headers for model moment file
cell_moments_header_path = user_arguments["cell_moments_header_path"]
# Specifies path to output directory
output_dir = user_arguments["output_dir"]
# Specifies output file name
output_file_name = user_arguments["output_file_name"]
# _____ File path declarations END _____
moments = np.loadtxt(moments_file_path,delimiter=",")




CMI_raw_data_folder_name = user_arguments["CeeMI_raw_data_dir_prefix"]
CMI_raw_data_folder_location = user_arguments["CeeMI_raw_data_dir_location"]
CMI_crm_data_folder_location = user_arguments["CeeMI_crm_data_dir_location"]
CMI_crm_data_folder_name = user_arguments["CeeMI_crm_data_dir_prefix"]


if not os.path.exists(CMI_raw_data_folder_location):
	os.mkdir(f"{CMI_raw_data_folder_location}")
if not os.path.exists(CMI_crm_data_folder_location):
	os.mkdir(f"{CMI_crm_data_folder_location}")
# parameter_dict_fi = {"foxo": foxo, "igfr": igfr}



cf.create_dir_of_cell_data_files({}, moments, f"{CMI_raw_data_folder_location}/{CMI_raw_data_folder_name}", cell_data_column_header_name_and_location=cell_moments_header_path, all_cells_are_conditioned_on=True)
try:
	os.mkdir(f"{CMI_crm_data_folder_location}/{CMI_crm_data_folder_name}")
except OSError as error:
	delete_file = input(f"{CMI_crm_data_folder_location}/{CMI_crm_data_folder_name} already exists, remove? (y/n)")
	if delete_file.lower() == "y":
		shutil.rmtree(f"{CMI_crm_data_folder_location}/{CMI_crm_data_folder_name}")
		os.mkdir(f"{CMI_crm_data_folder_location}/{CMI_crm_data_folder_name}")
	else:
		raise OSError
for cell_file_name in os.listdir(f"{CMI_raw_data_folder_location}/{CMI_raw_data_folder_name}"):
	prc.moment_file_to_crm_file(f"{CMI_raw_data_folder_location}/{CMI_raw_data_folder_name}/{cell_file_name}",f"{CMI_crm_data_folder_location}/{CMI_crm_data_folder_name}/{cell_file_name}", assumed_distribution="gamma", discretization_parameter=0.2)


cmi_all_cell, input_cc = mi.pcmi_at_cc(f"{CMI_crm_data_folder_location}/{CMI_crm_data_folder_name}",return_all_cell_values = True,sort_filename=1)

if not os.path.exists(output_dir):
	os.mkdir(f"{output_dir}")
np.savetxt(f"{output_dir}/single_cell_mi_array.csv", cmi_all_cell, delimiter=",")
np.savetxt(f"{output_dir}/cc_in.csv", input_cc, delimiter=",")
