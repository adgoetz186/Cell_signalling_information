import os
import numpy as np
import matplotlib.pyplot as plt

def population_moments_from_file(input_dose_folder_path,processed_data_filename,time_points,moments = 2):
	data_time = np.loadtxt(input_dose_folder_path+processed_data_filename,delimiter=",")
	time = data_time[0, :]
	time_indicies = np.searchsorted(time, np.array(time_points))
	data = data_time[1:, :]
	data = data[:,time_indicies]
	population_moments = np.zeros((moments,len(time_indicies)))
	for moment in range(1,moments+1):
		population_moments[moment-1] = np.average(data**moment,axis = 0)
	return population_moments


def population_moments_from_folder_of_dose_files(input_dose_folder_path,time_points,moments = 2):
	file_moment_list = []
	for dose_filename in os.listdir(input_dose_folder_path):
		print(dose_filename)
		file_moment_list.append(population_moments_from_file(input_dose_folder_path,dose_filename,time_points=time_points))
	moment_response_array = np.zeros((moments,len(time_points)*len(file_moment_list)))
	for moment_ind in range(moments):
		for file_ind in range(len(file_moment_list)):
			moment_response_array[moment_ind,file_ind*len(time_points):(file_ind+1)*len(time_points)] = file_moment_list[file_ind][moment_ind]
	header_title = []
	for moment in range(moments):
		for file_ind in range(len(file_moment_list)):
			for time in time_points:
				entry_dict = {}
				entry_dict["moment"] = f"{moment + 1}"
				entry_dict["dose"] = f"{os.listdir(input_dose_folder_path)[file_ind].split('.')[0]}"
				entry_dict["time"] = f"{time}_min"
				header_title.append(entry_dict)
	moment_response_array = np.reshape(moment_response_array,(1,-1))
	return moment_response_array, header_title


def moments_and_header_from_dose_folder(input_dose_folder_location,output_file_location,output_file_name,time_points = [],moments = 2):
	data,header_list = population_moments_from_folder_of_dose_files(input_dose_folder_location, time_points=time_points, moments=2)
	np.savetxt(output_file_location+output_file_name+".csv",data,delimiter=",")
	with open(output_file_location+output_file_name + "_header.txt","w") as header_file:
		for header in header_list:
			header_file.write(str(header)+"\n")

def single_cell_moments_from_steady_state_distributions(input_dose_folder_path,processed_data_filename,output_file_location,output_file_name,time_point_arrays,dose_name_list,moments = 2):
	response_traj = np.loadtxt(input_dose_folder_path + processed_data_filename + ".csv",delimiter=",")
	time_traj = response_traj[0,:]
	response_traj = response_traj[1:,:]
	cell_moment_array = np.zeros((np.shape(response_traj)[0], moments*len(time_point_arrays)))
	print(response_traj)
	for time_point_array_ind in range(len(time_point_arrays)):
		print(np.searchsorted(time_traj, time_point_arrays[time_point_array_ind]))
		cell_response_array = response_traj[:,np.searchsorted(time_traj, time_point_arrays[time_point_array_ind])]
		print(np.mean(cell_response_array,axis=1))
		print(np.var(cell_response_array,axis=1))
		for i in range(moments):
			cell_moment_array[:,time_point_array_ind+i*(len(time_point_arrays))] = np.mean(cell_response_array**(i+1),axis=1)
	np.savetxt(output_file_location+output_file_name+".csv",cell_moment_array,delimiter=",")
	
	

	with open(output_file_location+output_file_name + "_header.txt",'w') as header_file:
		for moment in range(moments):
			for dose_ind in range(len(dose_name_list)):
				entry_dict = {}
				entry_dict["moment"] = f"{moment + 1}"
				entry_dict["dose"] = f"{dose_name_list[dose_ind]}"
				header_file.write(str(entry_dict)+"\n")
	
