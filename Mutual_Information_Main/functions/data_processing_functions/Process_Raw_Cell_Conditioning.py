import numpy as np
import scipy.stats as st
import ast


def create_list_of_cell_data_files(conditioning_parameter_dict, cell_data_to_use, conditional_tolerance=0,
                                   all_cells_are_conditioned_on=False):
	"""Creates folder full of cells or cells with some aspect of cell state averaged over which can be used to obtain
	CeeMI, MI of CSAR, or something in between

	Input:

		-conditioning_parameter_dict (dict): dictionary with each key giving a single cell parameter (like initial igfr) and
		each value being a list of every cells value for that parameter, order is important

		-cell_data_to_use (ndarray): array with each row corresponding to a cell, the columns should be relevant values
		for defining the CRM for the cell (this would typically be moments or parameters)

		-conditional_tolerance (float): Cells with identical conditioning parameters will have their responses
		averaged over, this value allows cells which are almost identical to be combined.

		-all_cells_are_conditioned_on (boolean): A speed up for CMI, if listed as true the program will know to not
		aggregate any cells and so can jump right in to making files for every cell

	Output:

		-list_of_channels (list): A list with each entry contains parameters required to define a unique channel"""
	
	# If the below is true this handles special case which arises from no conditioning
	if len(conditioning_parameter_dict) == 0 and not all_cells_are_conditioned_on:
		if len(np.shape(cell_data_to_use)) == 2:
			cell_data_to_use = np.average(cell_data_to_use, axis=0)
		return [np.reshape(cell_data_to_use, (1, -1))]
	else:
		list_of_channels = []
		# conditioning_array: array which contains each cells set of values which are being conditioned on,
		# cells with identical conditioning values will be aggregated
		if not all_cells_are_conditioned_on:
			conditioning_array = np.zeros(
				(len(conditioning_parameter_dict), len(list(conditioning_parameter_dict.values())[0])))
			for conditional_parameter_index in range(len(conditioning_parameter_dict)):
				conditioning_array[conditional_parameter_index] = conditioning_parameter_dict[
					list(conditioning_parameter_dict.keys())[conditional_parameter_index]]
			conditioning_array = np.transpose(conditioning_array)
			cell_index_list = []
			# creates array where each cell is assigned an index corresponding to the first cell in the array with
			# identical conditioning parameters. If all cells are conditioned on, each cell is given a unique index
			if not all_cells_are_conditioned_on:
				for i in range(np.shape(conditioning_array)[0]):
					for j in range(i + 1):
						if np.all(conditioning_array[i] <= conditioning_array[j] + conditional_tolerance) and np.all(
								conditioning_array[i] >= conditioning_array[j] - conditional_tolerance):
							cell_index_list.append(j)
							break
		else:
			cell_index_list = range(np.shape(cell_data_to_use)[0])
		cell_index_list = np.array(cell_index_list)
		for i in range(np.max(cell_index_list) + 1):
			index_cell_list = list(np.where(cell_index_list == i)[0])
			
			cell_conditioning_values = {}
			for key in conditioning_parameter_dict.keys():
				cell_conditioning_values[key] = conditioning_parameter_dict[key][index_cell_list[0]]
			list_of_channels.append(cell_data_to_use[index_cell_list, :])
		return list_of_channels

def moment_array_to_crm(moment_array, assumed_distribution="gamma", discretization_parameter=0.05, percentile_cutoff=0.0005):
	"""Converts an array of moments into a conditional response matrix
			
			Input:
				-moment_array (ndarray): array of moments, columns index different doses and moments, should be ordered
				as follows, moment 1 dose 1, ..., moment 1 dose n, ..., moment 2 dose n.
				rows can be used to index cells, this code will obtain the average response over all provided cells
			
				-assumed_distribution (string): Specifies the assumed shape of the distribution of the conditional
				responses
				
				-discretization_parameter (float): Used to discretize continuous distributions using uniform binning,
				bin size is the inter quartile range of the smallest conditional response multiplied by the
				"discretization_parameter". It is approximately 1/n where n is how many bins will be given to the IQR
				of the narrowest conditional response. Uniform binning is essential as entries within a column of the
				conditional response matrix MUST refer to the same responses.
			
				-percentile_cutoff (float): specifies the effective min/max response by specifying the lowest percentile
				which should be included. Also removes the percentiles above "1-percentile_cutoff"
			
			Output:
				-conditional_response_matrix (ndarray): Let nS be the number of signals and nR be the number of
				responses the conditional_probability_matrix is an nS X nR array with element i,j specifying the
				probability of the jth response to the ith input"""
	lc_assumed_dist = assumed_distribution.lower()
	continuous_distribution_moment_dictionary = {"gamma":2,"gaussian":2}
	discrete_distribution_moment_dictionary = {"poisson":1,"nbinom":2}
	if lc_assumed_dist in continuous_distribution_moment_dictionary.keys():
		is_cont = True
	else:
		is_cont = False
	if is_cont:
		number_of_moments = continuous_distribution_moment_dictionary[lc_assumed_dist]
	else:
		number_of_moments = discrete_distribution_moment_dictionary[lc_assumed_dist]

	signal_count = int(np.shape(moment_array)[1] / number_of_moments)

	if len(np.shape(moment_array)) == 2:
		moments = np.average(moment_array, axis=0)
	else:
		moments = moment_array
	
	first_moments = moments[:signal_count]
	
	if number_of_moments > 1:
		second_moments = moments[signal_count:]
		variances = second_moments - first_moments ** 2

	
	min_list = []
	max_list = []
	iqr = []

	for signal_number in range(signal_count):
		mean = first_moments[signal_number]
		
		if lc_assumed_dist == "poisson":
			min_list.append(int(st.poisson.ppf(percentile_cutoff, mean)))
			max_list.append(int(st.poisson.ppf(1-percentile_cutoff, mean)))
			
		if lc_assumed_dist == "nbinom":
			variance = variances[signal_number]
			if variance < mean:
				raise ValueError(f"variance can not be bigger than mean for nbinom ({variance}>{mean})")
			p = mean/variance
			n = mean * p / (1 - p)
			min_list.append(int(st.nbinom.ppf(percentile_cutoff, n,p)))
			max_list.append(int(st.nbinom.ppf(1-percentile_cutoff, n,p)))
		
		if lc_assumed_dist == "gamma":
			variance = variances[signal_number]
			scale0 = variance / mean
			shape0 = mean / scale0
			min_list.append(st.gamma.ppf(percentile_cutoff, shape0, scale=scale0))
			max_list.append(st.gamma.ppf(1-percentile_cutoff, shape0, scale=scale0))
			iqr.append(st.gamma.ppf(.75, shape0, scale=scale0) - st.gamma.ppf(.25, shape0, scale=scale0))
			
		if lc_assumed_dist == "gaussian":
			variance = variances[signal_number]
			min_list.append(st.norm.ppf(percentile_cutoff, loc = mean, scale=variance**.5))
			max_list.append(st.norm.ppf(1-percentile_cutoff,  loc = mean, scale=variance**.5))
			iqr.append(st.norm.ppf(.75,  loc = mean, scale=variance**.5) - st.norm.ppf(.25,  loc = mean, scale=variance**.5))
		# Additional continuous distribution forms can go here
		
	# Sets min and max allowed receptor to be such that at worst only
	# "2*percentile_cutoff * 100%" of any response will be neglected
	min_value = min(min_list)
	max_value = max(max_list)
	
	if is_cont:
		# Bin size for discretization, based on the smallest response's interquartile range
		# For a Discretization_Parameter of 0.2, 50% of the sharpest peaked response gets 5 bins
		bin_size = min(iqr) * discretization_parameter
		discretized_receptor_counts = np.arange(min_value, max_value + bin_size, bin_size)
		CRM = np.zeros((signal_count, np.size(discretized_receptor_counts)-1))
	else:
		print(max_value,min_value)
		CRM = np.zeros((signal_count, max_value-min_value+1))

	for signal_number in range(signal_count):
		mean = first_moments[signal_number]
		if lc_assumed_dist == "poisson":
			values_to_use = np.arange(min_list[signal_number],max_list[signal_number]+1)
			CRM[signal_number] = np.concatenate((np.zeros(min_list[signal_number] - min_value),st.poisson.pmf(values_to_use, mean),np.zeros(max_value - max_list[signal_number])))
		
		if lc_assumed_dist == "nbinom":
			variance = variances[signal_number]
			p = mean / variance
			n = mean * p / (1 - p)
			values_to_use = np.arange(min_list[signal_number],max_list[signal_number]+1)
			CRM[signal_number] = np.concatenate((np.zeros(min_list[signal_number] - min_value),st.nbinom.pmf(values_to_use, n,p),np.zeros(max_value - max_list[signal_number])))
		
		if lc_assumed_dist == "gamma":
			variance = variances[signal_number]
			scale0 = variance / mean
			shape0 = mean / scale0
			# Obtains the integral of the probability within the bin edges
			cdf_values = st.gamma.cdf(discretized_receptor_counts, shape0, scale=scale0)
			CRM[signal_number] = cdf_values[1:]-cdf_values[:-1]
			# Normalize the CRM (prior to normalization sum of distribution ~1-2*percentile_cutoff)
			CRM[signal_number] /= np.sum(CRM[signal_number])
		if lc_assumed_dist == "gaussian":
			variance = variances[signal_number]
			cdf_values = st.norm.cdf(discretized_receptor_counts, loc = mean, scale=variance**.5)
			CRM[signal_number] = cdf_values[1:]-cdf_values[:-1]
			CRM[signal_number] /= np.sum(CRM[signal_number])
		# Additional continuous distribution forms can go here
	return CRM

def moment_list_to_crm_list(moment_list, assumed_distribution ="gamma", discretization_parameter=0.05, percentile_cutoff=0.0005):
	crm_list = []
	for i in moment_list:
		crm_list.append(moment_array_to_crm(i,assumed_distribution = assumed_distribution, discretization_parameter = discretization_parameter, percentile_cutoff=percentile_cutoff))
	return crm_list


