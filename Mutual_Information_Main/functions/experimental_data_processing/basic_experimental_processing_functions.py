import numpy as np
from matplotlib import pyplot as plt
import copy as cp
import os


def load_trajectories(file_name):
	full_array = np.loadtxt(file_name, delimiter=",")
	Time_Values = full_array[:, 0]
	main_array = np.transpose(full_array[:, 1:])
	return Time_Values, main_array


def remove_background(trajectory_array, background_file_name):
	Background_Array = np.loadtxt(background_file_name, delimiter=",")
	trajectory_array = np.transpose(trajectory_array)
	Background = Background_Array[:, 1]
	# Removes the background and returns the trajectory array to having time on the x axis
	trajectory_array = np.transpose(
		trajectory_array - np.transpose(np.tile(Background, (np.shape(trajectory_array)[1], 1))))
	return trajectory_array


def remove_timepoints(trajectory_array, time_array, array_of_time_indexes_to_remove):
	time_indexes_to_keep = list(
		set([i for i in range(np.shape(trajectory_array)[1])]) - set(array_of_time_indexes_to_remove))
	return time_array[time_indexes_to_keep], trajectory_array[:, time_indexes_to_keep]


def remove_offset(list_of_trajectory_arrays):
	avg_signal_t_0 = 0
	total_divide_value = 0
	for i in list_of_trajectory_arrays:
		avg_signal_t_0 += np.sum(i[:, 0])
		total_divide_value += np.size(i[:, 0])
	avg_signal_t_0 /= total_divide_value
	list_of_shifts = []
	for i in list_of_trajectory_arrays:
		list_of_shifts.append(avg_signal_t_0 - np.average(i[:, 0]))
	shifted_list_of_Trajectory_arrays = []
	for i in range(len(list_of_trajectory_arrays)):
		shifted_list_of_Trajectory_arrays.append(list_of_trajectory_arrays[i] + list_of_shifts[i])
	return shifted_list_of_Trajectory_arrays


def relative_response(trajectory_array):
	# run this after any main data processing
	for i in range(np.shape(trajectory_array)[0]):
		trajectory_array[i] /= trajectory_array[i, 0]
	return trajectory_array


def reweight_list_of_trajectories(trajectory_list, additional_constant_mult=1.0):
	t0val = 0
	tdiv = 0
	for i in trajectory_list:
		t0val += np.sum(i[:, 0])
		tdiv += np.size(i[:, 0])
	t0val /= tdiv
	for i in range(len(trajectory_list)):
		trajectory_list[i] *= (additional_constant_mult / t0val)
	return trajectory_list

