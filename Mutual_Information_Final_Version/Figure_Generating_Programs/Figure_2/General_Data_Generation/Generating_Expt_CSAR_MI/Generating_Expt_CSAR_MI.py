import numpy as np
import os
import Mutual_Information_Final_Version.functions.data_processing_functions.Conditioning_Functions as cf
import Mutual_Information_Final_Version.functions.data_processing_functions.Process_Raw_Cell_Conditioning as prc
import Mutual_Information_Final_Version.functions.mutual_information_functions.MI_Calculation_Functions as mi
import Mutual_Information_Final_Version.functions.file_processing_functions.user_input_functions as uif
import Mutual_Information_Final_Version.functions.file_processing_functions.make_folder as mk

# Generates mutual information value of cell state agnostic response from experimental data for EGFR

# _____ Load arguments BEGIN _____
# Specifies path to file which specifies user specified arguments
user_input_file_path = "/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Figure_Generating_Programs/Figure_2/General_Data_Generation/Generating_Expt_CSAR_MI/args.txt"
user_arguments = uif.read_file_to_arg_dict(user_input_file_path)
# _____ Load arguments END _____


# _____ File path declarations BEGIN _____
# Specifies path to experimental moment files
expt_mom_1 = user_arguments["expt_mom_1"]
expt_mom_2 = user_arguments["expt_mom_2"]

# Output file location
output_information_dir = user_arguments["output_information_dir"]

# Intermediate folder locations/names
CeeMI_raw_data_folder_location = user_arguments["CeeMI_raw_data_folder_location"]
CeeMI_crm_data_folder_location = user_arguments["CeeMI_crm_data_folder_location"]
CeeMI_raw_data_folder = user_arguments["CeeMI_raw_data_folder"]
CeeMI_crm_data_folder = user_arguments["CeeMI_crm_data_folder"]
# _____ File path declarations END _____


# _____ Loading files BEGIN _____
# loads moment files
expt_mom_1_array = np.load(expt_mom_1)
expt_mom_2_array = np.load(expt_mom_2)
# _____ Loading files End _____

# _____ Main code BEGIN _____
# conversion factor from a.u. to receptor count
au_cf = user_arguments["au_cf"]

# converts moments from a.u. into appropriate units
expt_mom_1_array/=au_cf
expt_mom_2_array/=au_cf**2

# combines moment arrays
moment_array = np.hstack((expt_mom_1_array,expt_mom_2_array))

# creates raw cell conditioning folder and fills it with raw cell files (moments)
cf.create_dir_of_cell_data_files({}, moment_array,f"{CeeMI_raw_data_folder_location}/{CeeMI_raw_data_folder}")

# creates and populates folder of conditional response matrices
mk.make_folder(f"{CeeMI_crm_data_folder_location}/{CeeMI_crm_data_folder}")
for cell_file_name in os.listdir(f"{CeeMI_raw_data_folder_location}/{CeeMI_raw_data_folder}"):
	prc.moment_file_to_crm_file(f"{CeeMI_raw_data_folder_location}/{CeeMI_raw_data_folder}/{cell_file_name}",f"{CeeMI_crm_data_folder_location}/{CeeMI_crm_data_folder}/{cell_file_name}", assumed_distribution=user_arguments["distribution_type"],discretization_parameter=user_arguments["discretization_parameter"])

# obtains channel capacity for MI of CSAR
cmi_all_cell, cmi_cc_in = mi.pcmi_at_cc(f"{CeeMI_crm_data_folder_location}/{CeeMI_crm_data_folder}", return_all_cell_values=True)

# saves output
if not os.path.exists(output_information_dir):
	os.mkdir(output_information_dir)
np.savetxt(f"{output_information_dir}/single_cell_MI.csv",cmi_all_cell,delimiter=",")
np.savetxt(f"{output_information_dir}/cc_input_distribution",cmi_cc_in,delimiter=",")
# _____ Main code End _____
