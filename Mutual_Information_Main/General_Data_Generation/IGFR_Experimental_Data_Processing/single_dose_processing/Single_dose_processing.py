import os
import numpy as np
import Mutual_Information_Final_Version.functions.experimental_data_processing.basic_experimental_processing_functions as epf
import Mutual_Information_Final_Version.functions.file_processing_functions.user_input_functions as uif

# Processes the raw single dose florescence datasets provided by Heiser Lab
# https://pubmed.ncbi.nlm.nih.gov/31838146/

# _____ Load arguments BEGIN _____
# Specifies path to file which specifies user specified arguments
user_input_file_path = "no_250_dose_args.txt"
user_arguments = uif.read_file_to_arg_dict(user_input_file_path)
# _____ Load arguments END _____

# _____ File path declarations BEGIN _____
# location of experimental data folder as provided by the Heiser lab
raw_dose_folder = user_arguments["raw_dose_folder"]
# location of background fluorescence data file as provided by the Heiser lab
background_file = user_arguments["background_file"]
# location of folder for processed experimental data
output_dose_folder = user_arguments["output_dose_folder"]
# _____ File path declarations END _____

# _____ Loading files BEGIN _____
list_of_inputs = user_arguments["list_of_inputs"]
trajectory_list = []
for i in list_of_inputs:
	times, traj = epf.load_trajectories(f"{raw_dose_folder}{i}.csv")
	trajectory_list.append(traj)
# _____ Loading files END _____

# _____ Main code BEGIN _____
for i in range(len(trajectory_list)):
	# Removes background florescence
	trajectory_list[i] = epf.remove_background(trajectory_list[i],background_file)
	
	# Removes the data collected at 3 minutes, which is corrupted
	new_times, trajectory_list[i] = epf.remove_timepoints(trajectory_list[i], times, [1])
	
# Shifts each measured response so the initial t=0 values agree with one another
trajectory_list = epf.remove_offset(trajectory_list)

# Converts values from a.u. to protein count
trajectory_list = epf.reweight_list_of_trajectories(trajectory_list, additional_constant_mult=user_arguments["foxo_0"])

# Saves responses
for i in range(len(trajectory_list)):
	experimental_trajectory = np.vstack((new_times, trajectory_list[i]))
	np.savetxt(f"{output_dose_folder}{list_of_inputs[i]}", experimental_trajectory, delimiter=",")
# _____ Main code END _____
