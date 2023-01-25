import numpy as np
import ast
import Mutual_Information_Final_Version.functions.file_processing_functions.make_folder as mk

def create_dir_of_cell_data_files(conditioning_parameter_dict, cell_data_to_save, dir_path, cell_data_column_header_name_and_location="", conditional_tolerance=0,all_cells_are_conditioned_on=False):
	"""Creates folder full of cells or cells with some aspect of cell state averaged over which can be used to obtain
	CeeMI, MI of CSAR, or something in between
	
	Input:
	
		-conditioning_parameter_dict (dict): dictionary with each key giving a single cell parameter (like initial igfr) and
		each value being a list of ever cells value for that parameter, order is important
		
		-cell_data_to_save (ndarray): array with each row corresponding to a cell, the columns should be relevant values
		for defining the CRM for the cell (this would typically be moments or parameters)
		
		-dir_path (str): The path and name of the directory
		
		-cell_data_column_header_name_and_location (list): list of dictionaries where each dictionary specifies
		identifying information for the corresponding column of the cell_data_to_save
		
		-conditional_tolerance (float): Cells with identical conditioning parameters will have their responses
		averaged over, this value allows cells which are almost identical to be combined.
		
		-all_cells_are_conditioned_on (boolean): A speed up for CMI, if listed as true the program will know to not
		aggregate any cells and so can jump right in to making files for every cell
		
	Output:
	
		-returns none but generates a folder of the specified name and at the specified location
		folder contains a unique file for each unique value of parameters in the conditioning_parameter_dict
		each file contains an array which is a subset of cell_data_to_save,
		with each file having the full data of some number of cells
		it should be possible to obtain the population level conditional response matrix for the cells in each bin"""
	# Obtains column specifying data (input dose, time of dose, etc) from file if file is specified
	if cell_data_column_header_name_and_location != "":
		with open(cell_data_column_header_name_and_location) as cell_data_column_header_file:
			cell_data_column_header_list = [ast.literal_eval(i) for i in cell_data_column_header_file.readlines()]
		cell_data_column_header = cell_data_column_header_list
	else:
		cell_data_column_header = ""
	
	# creates a folder at the specified path
	mk.make_folder(dir_path)
	# If the below is true this handles special case which arises from no conditioning
	if len(conditioning_parameter_dict) == 0 and not all_cells_are_conditioned_on:
		if cell_data_column_header == "":
			np.savetxt(f"{dir_path}/population_response", cell_data_to_save, delimiter=",")
		else:
			np.savetxt(f"{dir_path}/population_response", cell_data_to_save, delimiter=",", header=f"{cell_data_column_header}")
		return
	else:
		# conditioning_array: array which contains each cells set of values which are being conditioned on,
		# cells with identical conditioning values will be aggregated
		if not all_cells_are_conditioned_on:
			conditioning_array = np.zeros((len(conditioning_parameter_dict),len(list(conditioning_parameter_dict.values())[0])))
			for conditional_parameter_index in range(len(conditioning_parameter_dict)):
				conditioning_array[conditional_parameter_index] = conditioning_parameter_dict[list(conditioning_parameter_dict.keys())[conditional_parameter_index]]
			conditioning_array = np.transpose(conditioning_array)
			cell_index_list = []
			# creates array where each cell is assigned an index corresponding to the first cell in the array with
			# identical conditioning parameters. If all cells are conditioned on, each cell is given a unique index
			if not all_cells_are_conditioned_on:
				for i in range(np.shape(conditioning_array)[0]):
					for j in range(i+1):
						if np.all(conditioning_array[i]<=conditioning_array[j]+conditional_tolerance) and np.all(conditioning_array[i]>=conditioning_array[j]-conditional_tolerance):
							cell_index_list.append(j)
							break
		else:
			cell_index_list = range(np.shape(cell_data_to_save)[0])
		cell_index_list = np.array(cell_index_list)
		for i in range(np.max(cell_index_list)+1):
			index_cell_list = list(np.where(cell_index_list == i)[0])
			if len(index_cell_list) == 1:
				cell_name = f"cell_{i}"
			else:
				cell_name = f"averaged_cell_{i}"
			cell_conditioning_values = {}
			for key in conditioning_parameter_dict.keys():
				cell_conditioning_values[key] = conditioning_parameter_dict[key][index_cell_list[0]]
			if cell_data_column_header == "":
				np.savetxt(f"{dir_path}/{cell_name}", cell_data_to_save[index_cell_list, :], delimiter=",",
				           header=f"{cell_conditioning_values}")
			else:
				np.savetxt(f"{dir_path}/{cell_name}", cell_data_to_save[index_cell_list,:], delimiter=",", header=f"{cell_conditioning_values}\n{cell_data_column_header}")
		return