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
user_input_file_path = "/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Figure_Generating_Programs/Figure_3/General_Data_Generation/Generating_model_CeeMI_trajectory/args.txt"
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
all_parameter_list = list(np.loadtxt(cell_params_header_path,dtype=str))
cell_index_list = list(range(np.shape(moments)[0]))

number_of_cells = user_arguments["number_of_cells"]
number_of_runs = user_arguments["number_of_runs"]
number_of_inputs = user_arguments["number_of_inputs"]
time_list = user_arguments["time_list"]
single_cell_mi_array = np.zeros((number_of_runs*number_of_cells,len(time_list)))

ceemi_array = np.zeros((number_of_runs,len(time_list)))
cc_ceemi_array = np.zeros((number_of_runs,number_of_inputs))
for time_index in range(len(time_list)):
    cell_count = 0
    for run_num in range(number_of_runs):
        random.shuffle(cell_index_list)
        selected_cells = cell_index_list[:number_of_cells]

        CMI_raw_data_folder_name = f"{user_arguments['CeeMI_raw_data_dir_prefix']}_run_{run_num}_time_{time_list[time_index]}"
        CMI_raw_data_folder_location = user_arguments['CeeMI_raw_data_dir_location']
        CMI_crm_data_folder_location = user_arguments['CeeMI_crm_data_dir_location']
        CMI_crm_data_folder_name = f"{user_arguments['CeeMI_crm_data_dir_prefix']}_run_{run_num}_time_{time_list[time_index]}"


        params_select = params[selected_cells,:]
        
        moments_select = moments[selected_cells,:]
        list_of_param_values = [10 ** params_select[:, i] for i in range(np.shape(params_select)[1])]
        param_dict = dict(zip(all_parameter_list, list_of_param_values))
        
        if not os.path.exists(CMI_raw_data_folder_location):
            os.mkdir(f"{CMI_raw_data_folder_location}")
        if not os.path.exists(CMI_crm_data_folder_location):
            os.mkdir(f"{CMI_crm_data_folder_location}")

        cf.create_dir_of_cell_data_files(param_dict, moments_select, f"{CMI_raw_data_folder_location}/{CMI_raw_data_folder_name}",
                                            cell_data_column_header_name_and_location=cell_moments_header_path,all_cells_are_conditioned_on=True)
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
        for cell_file_name in os.listdir(f"{CMI_raw_data_folder_location}/{CMI_raw_data_folder_name}"):
            prc.moment_file_to_crm_file(f"{CMI_raw_data_folder_location}/{CMI_raw_data_folder_name}/{cell_file_name}",f"{CMI_crm_data_folder_location}/{CMI_crm_data_folder_name}/{cell_file_name}",
                                        lock_parameter={"time": f"{time_list[time_index]}_min"},
                                        assumed_distribution=user_arguments["distribution_type"], discretization_parameter=user_arguments["discretization_parameter"])
        cmi_all_cell, cmi_cc_in = mi.pcmi_at_cc(f"{CMI_crm_data_folder_location}/{CMI_crm_data_folder_name}",return_all_cell_values=True)
        ceemi_array[run_num,time_index] = np.average(cmi_all_cell)
        for i in range(np.shape(cmi_all_cell)[0]):
            single_cell_mi_array[cell_count,time_index] = cmi_all_cell[i]
            cell_count += 1
        if time_list[time_index] == user_arguments["time_for_cc_input"]:
            cc_ceemi_array[run_num] = cmi_cc_in
if not os.path.exists(output_dir):
    os.mkdir(f"{output_dir}")
np.savetxt(f"{output_dir}/cc_cmi_array.csv",cc_ceemi_array,delimiter=",")
np.savetxt(f"{output_dir}/cmi_array.csv",ceemi_array,delimiter=",")
np.savetxt(f"{output_dir}/single_cell_mi_array.csv",single_cell_mi_array,delimiter=",")
# _____ Main code END _____
