import random
import numpy as np
import os
import shutil
import Mutual_Information_Final_Version.functions.data_processing_functions.Conditioning_Functions as cf
import Mutual_Information_Final_Version.functions.data_processing_functions.Process_Raw_Cell_Conditioning as prc
import Mutual_Information_Final_Version.functions.mutual_information_functions.MI_Calculation_Functions as mi
import Mutual_Information_Final_Version.functions.file_processing_functions.user_input_functions as uif

# This program generates the model CeeMI trajectories for the igfr/akt/foxo pathway

# _____ Load arguments BEGIN _____
# Specifies path to file which specifies user specified arguments
user_input_file_path = "/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Figure_Generating_Programs/Figure_3/General_Data_Generation/Generating_Model_Step_Dose_CeeMI/args.txt"
user_arguments = uif.read_file_to_arg_dict(user_input_file_path)
# _____ Load arguments END _____

# _____ File path declarations BEGIN _____
# Specifies path to drawn model cell parameters
params_file_path = user_arguments["params_file_path"]
# Specifies path to drawn model cell moments
moments_file_path = user_arguments["moments_file_path"]
# Specifies path to headers for model parameter file
cell_params_header_path = user_arguments["cell_params_header_path"]
# Specifies path to headers for model moment file
cell_moments_header_path = user_arguments["cell_moments_header_path"]
# Specifies path to output directory
output_dir = user_arguments["output_dir"]
# _____ File path declarations END _____

# _____ Loading files BEGIN _____
# Loads parameters of drawn cells
params = np.loadtxt(params_file_path,delimiter=",")

# Loads moments of drawn cells
moments = np.loadtxt(moments_file_path,delimiter=",")
# _____ Loading files END _____


# _____ Main code BEGIN _____
even_tpl = tuple([2*i for i in range(np.shape(moments)[1]//2)])
odd_tpl = tuple([2*i+1 for i in range(np.shape(moments)[1]//2)])

# Averages 60 and 90 minute moments
moments = (moments[:,even_tpl] + moments[:,odd_tpl])/2

all_parameter_list = list(np.loadtxt(cell_params_header_path,dtype=str))
cell_index_list = list(range(np.shape(moments)[0]))

number_of_cells = user_arguments["number_of_cells"]
number_of_inputs = user_arguments["number_of_inputs"]
random.shuffle(cell_index_list)
selected_cells = cell_index_list[:number_of_cells]
params_select = params[selected_cells,:]
single_cell_mi_array = np.zeros_like(moments)
single_cell_mi_params = np.zeros_like(params)

cc_ceemi_array = np.zeros((1,number_of_inputs))


CMI_raw_data_folder_name = user_arguments['CeeMI_raw_data_dir_prefix']
CMI_raw_data_folder_location = user_arguments['CeeMI_raw_data_dir_location']
CMI_crm_data_folder_location = user_arguments['CeeMI_crm_data_dir_location']
CMI_crm_data_folder_name = user_arguments['CeeMI_crm_data_dir_prefix']



moments_select = moments[selected_cells,:]
print(np.shape(moments_select))
list_of_param_values = [10 ** params_select[:, i] for i in range(np.shape(params_select)[1])]
param_dict = dict(zip(all_parameter_list, list_of_param_values))

if not os.path.exists(CMI_raw_data_folder_location):
    os.mkdir(f"{CMI_raw_data_folder_location}")
if not os.path.exists(CMI_crm_data_folder_location):
    os.mkdir(f"{CMI_crm_data_folder_location}")

cf.create_dir_of_cell_data_files(param_dict, moments_select, f"{CMI_raw_data_folder_location}/{CMI_raw_data_folder_name}",all_cells_are_conditioned_on=True)
try:
    os.mkdir(f"{CMI_crm_data_folder_location}/{CMI_crm_data_folder_name}")
except OSError as error:
    delete_file = input(
        f"{CMI_crm_data_folder_location}/{CMI_crm_data_folder_name} already exists, remove? (y/n)")
    if delete_file.lower() == "y":
        shutil.rmtree(f"{CMI_crm_data_folder_location}/{CMI_crm_data_folder_name}")
        os.mkdir(f"{CMI_crm_data_folder_location}/{CMI_crm_data_folder_name}")
    else:
        raise OSError
print(user_arguments["distribution_type"])
print(user_arguments["discretization_parameter"])
for cell_file_name in os.listdir(f"{CMI_raw_data_folder_location}/{CMI_raw_data_folder_name}"):
    prc.moment_file_to_crm_file(f"{CMI_raw_data_folder_location}/{CMI_raw_data_folder_name}/{cell_file_name}",f"{CMI_crm_data_folder_location}/{CMI_crm_data_folder_name}/{cell_file_name}",assumed_distribution=user_arguments["distribution_type"], discretization_parameter=user_arguments["discretization_parameter"])
cmi_all_cell,cc_in = mi.pcmi_at_cc(f"{CMI_crm_data_folder_location}/{CMI_crm_data_folder_name}",return_all_cell_values = True,sort_filename=1)

        
if not os.path.exists(output_dir):
    os.mkdir(f"{output_dir}")
np.savetxt(f"{output_dir}/mdl_moments_used.csv", moments_select, delimiter=",")
np.savetxt(f"{output_dir}/mdl_params_used.csv", params_select, delimiter=",")
np.savetxt(f"{output_dir}/cc_in.csv",cc_in,delimiter=",")
np.savetxt(f"{output_dir}/single_cell_mi_array.csv",cmi_all_cell,delimiter=",")
# _____ Main code END _____
