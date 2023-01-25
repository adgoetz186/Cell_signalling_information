import os
import numpy as np
import copy as cp
import Mutual_Information_Final_Version.functions.file_processing_functions.user_input_functions as uif

# Processes the raw multi-dose florescence data provided by Heiser Lab
# https://pubmed.ncbi.nlm.nih.gov/31838146/

# _____ Load arguments BEGIN _____
# Specifies path to file which specifies user specified arguments
user_input_file_path = "/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Experimental_Data_Processing/multidose_processing/args.txt"
user_arguments = uif.read_file_to_arg_dict(user_input_file_path)
# _____ Load arguments END _____

# _____ File path declarations BEGIN _____
# location for multidose file as provided by Heiser Lab
raw_mutlidose_file = user_arguments["raw_mutlidose_file"]
# location to store processed multidose file
output_mutlidose_file = user_arguments["output_mutlidose_file"]
# _____ File path declarations END _____

# _____ Loading files BEGIN _____
mtd_data = np.loadtxt(raw_mutlidose_file,delimiter=",")
# _____ Loading files End _____

# _____ Main code BEGIN _____
# Background measurements for each cell to remove
background = mtd_data[:,0]

# Time Values
times = mtd_data[:,1]

# Measured flourescence trajectories
trajectories = mtd_data[:,2:]

# Removes measured background values
background_removed_trajectories = trajectories - np.reshape(background,(-1,1))

# converts from a.u. to average protein count
background_removed_trajectories *= user_arguments["foxo_0"]/np.average(background_removed_trajectories[0])

# rearranges and saves array
background_removed_trajectories = np.transpose(background_removed_trajectories)
background_removed_trajectories = np.vstack((times,background_removed_trajectories))
np.savetxt(output_mutlidose_file,background_removed_trajectories,delimiter=',')
# _____ Main code END _____