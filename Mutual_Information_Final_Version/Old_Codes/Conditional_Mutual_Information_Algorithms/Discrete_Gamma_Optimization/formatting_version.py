import os
import numpy as np
from functools import partial
import time
import ast
from scipy.optimize import Bounds
from scipy.optimize import LinearConstraint
from scipy.optimize import minimize
import multiprocessing as mp
from Mutual_Information_Final_Version.Old_Codes.MI_Functions.Input_Functions import mutual_information_of_cell
import Mutual_Information_Final_Version.functions.data_processing_functions.Process_Raw_Cell_Conditioning as prc


# This method requires specifying a datafile of the following format:
#
# Datafile should have the moments of all inputs organized in the following manner:
# Each cell gets its own row. each moment will be given by uij where i designates the moment number (first moment, second moment)
# and j designates the signal input number (signal 1, signal 2,...)
# The following is for a system with m moments and n signals
# u11,u12,u13,...,u1n,u21,u22,u23,u2n,...,umn
#
# It may be good to make a program to directly generate moments from parameters, but care should be done to ensure bad
# signal values are not used (signals for which the assumed distribution is incorrect)

def mutual_information_for_all_cells(input_vector, data_list, averaged_output,step_size_multiplier = 0.25):
	# Takes in a 1D numpy array for an input vector as well as the data_list, the list containing all cell data
	# Handles the multiprocessing order
	# len(data_list)
	timestart = time.time()
	full_order_list = [i for i in range(len(data_list))]
	# full_order_list = [48]
	pool = mp.Pool(processes=mp.cpu_count())
	func = partial(mutual_information_of_cell, "gamma", input_vector, data_list,step_size_multiplier)
	ResB_P = pool.map(func, full_order_list)
	pool.close()
	if averaged_output:
		return np.average(ResB_P)
	else:
		return ResB_P



# print(directory_location)

def generate_reduced_moment_array(raw_file_path, raw_file_name, output_file_name="", lock_parameter={}, number_of_moments=2,
                            output_file_path=""):
	# This is going to need sorting
	# Currently this assumes the order of columns goes low dose -> high dose, low moment -> high moment
	# It would be good to remove this restriction, but not essential
	# It might also be useful to make the moment identifier a number rather than "first" or "second"
	# Input:
	#   lock_parameter: specifies a dictionary of keys and values which must be in a column's dictionary for the moment to be used
	#       Often this will be a single entry dictionary of the form "{time:i_min}"
	if raw_file_path[-1] == "/":
		raw_file_path = raw_file_path[:-1]
	with open(f"{raw_file_path}/{raw_file_name}", "r") as file:
		file_contents = file.readlines()
	header_string = ""
	
	if file_contents[1][0] != "#":
		column_index_list = ast.literal_eval(file_contents[0][2:])
	else:
		header_string += file_contents[0][2:-1]
		column_index_list = ast.literal_eval(file_contents[1][2:])
	moment_array = np.loadtxt(
		f"{raw_file_path}/{raw_file_name}", delimiter=",")
	column_indices_to_keep = []
	new_column_index_list = []
	for col_index in range(len(column_index_list)):
		add_entry = True
		for lp in lock_parameter.keys():
			if lock_parameter[lp] != column_index_list[col_index][lp]:
				add_entry = False
		if add_entry:
			column_indices_to_keep.append(col_index)
			new_column_index_list.append(column_index_list[col_index])
	column_index_list = new_column_index_list
	header_string += f"\n{column_index_list}"
	moment_array = moment_array[:, column_indices_to_keep]
	return
	


if __name__ == '__main__':
	os.chdir("../../../..")
	input_directory_location = "Mutual_Information_Final_Version/Data/IGFR/Raw_Cell_Files_Conditioning/MI_bn_5_cn_500/"
	#output_directory_location = "Data/EGFR/CC_Cond_MI/"
	
	for input_file_name in os.listdir(input_directory_location):
		# See top of script for instructions on data file format
		with open(input_directory_location + input_file_name) as filename:
			data_list = filename.readlines()
		print(data_list[0])
		# comments matter
		moment_count = 2
		# provide the number of moments used in the distribution you will use
		# right now this value might as well be locked in at 2 as only negative binomial fits are used
		# my hope is to extend this to gamma, poisson, and gaussian
		moment_array,header_data = prc.moment_file_to_moments(input_directory_location,input_file_name,lock_parameter={"time":"90_min"})

		#print(moment_array)
		signal_count = int(np.shape(moment_array)[1] / moment_count)
		bounds = Bounds(np.zeros(signal_count), np.ones(signal_count))
		linear_constraint = LinearConstraint(np.ones(signal_count), [1], [1])
		x0 = np.ones(signal_count) / signal_count
		output_file_name = input_file_name.split(".")[0] + ".txt"
		res = minimize(mutual_information_for_all_cells, x0, method='trust-constr', constraints=linear_constraint,
		               options={'verbose': 2}, bounds=bounds, args=(moment_array, True))
		output = np.reshape(res['x'], (1, -1))
		print(res)
		# input(mutual_information_from_matrix(x0,crm))
		
		#np.savetxt(f"{output_directory_location}CMI/{input_file_name.split('.')[0]}_CMI.csv", output, delimiter=",",
		#           header=str(res['fun'] * -1))
		#np.savetxt(f"{output_directory_location}sc_CMI/{input_file_name.split('.')[0]}_sc_CMI.csv",
		#           mutual_information_for_all_cells(res['x'], data_list, False))