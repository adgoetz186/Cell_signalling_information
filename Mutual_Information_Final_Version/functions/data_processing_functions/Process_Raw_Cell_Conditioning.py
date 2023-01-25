import numpy as np
import scipy.stats as st
import ast


def moment_array_to_crm(moment_array, assumed_distribution="gamma", discretization_parameter=5, percentile_cutoff=0.0005):
	"""Converts an array of moments into a conditional response matrix
			
			Input:
				-moment_array (ndarray): array of moments, columns index different doses and moments, should be ordered
				as follows, moment 1 dose 1, ..., moment 1 dose n, ..., moment 2 dose n.
				rows can be used to index cells, this code will obtain the average response over all provided cells
			
				-assumed_distribution (string): Specifies the assumed shape of the distribution of the conditional
				responses
				
				-discretization_parameter (int): Used to discretize continuous distributions using uniform binning,
				bin size is the inter quartile range of the smallest conditional response divided by the
				"discretization_parameter". It is approximately how many bins will be given to the IQR of the
				narrowest conditional response. Uniform binning is essential as entries within a column of the
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
	print(moments)
	
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


def moment_file_to_moments(file_path, lock_parameter=None):
	"""Converts a file containing moments into an array of moments
	
		Input:
		
			-file_path (str): path to file containing moments
			
			-lock_parameter (dict): specifies a dictionary of keys and values which must be in a column's dictionary
			for the moment to be used, often this will be a single entry dictionary of the form "{time:i_min}"
			
		Output:
		
			-array_of_moments (ndarray): array of moments with rows indexing unique cells and columns being described
			by the header string
			
			-header_string (string): describes the nature of each column in the "array_of_moments" """
	if lock_parameter is None:
		lock_parameter = {}
	with open(file_path, "r") as file:
		file_contents = file.readlines()
	header_string = ""
	
	if file_contents[1][0] != "#":
		column_index_list = ast.literal_eval(file_contents[0][2:])
	else:
		header_string += file_contents[0][2:-1]
		column_index_list = ast.literal_eval(file_contents[1][2:])
	moment_array = np.loadtxt(file_path, delimiter=",")
	if len(np.shape(moment_array)) == 1:
		moment_array = np.reshape(moment_array,(1,-1))
	if lock_parameter != {}:
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
	else:
		column_indices_to_keep = [i for i in range(np.shape(moment_array)[1])]
	
	array_of_moments = moment_array[:, column_indices_to_keep]
	return array_of_moments, header_string


def moment_file_to_crm_file(file_path, output_file_path, lock_parameter=None, assumed_distribution ="gamma", discretization_parameter=5, percentile_cutoff=0.0005):
	"""Converts a file containing moments into a file containing a conditional response matrix
	
	Input:
	
		-moment_array (ndarray): array of moments, columns index different doses and moments, should be ordered
		as follows, moment 1 dose 1, ..., moment 1 dose n, ..., moment 2 dose n.
		rows can be used to index cells, this code will obtain the average response over all provided cells
		
		-assumed_distribution (string): Specifies the assumed shape of the distribution of the conditional responses
		
		-discretization_parameter (int): Used to discretize continuous distributions using uniform binning,
		bin size is the inter quartile range of the smallest conditional response divided by the
		"discretization_parameter". It is approximately how many bins will be given to the IQR of the
		narrowest conditional response. Uniform binning is essential as entries within a column of the
		conditional response matrix MUST refer to the same responses.
		
		-percentile_cutoff (float): specifies the effective min/max response by specifying the lowest percentile
		which should be included. Also removes the percentiles above "1-percentile_cutoff"
		
	Output:
	
		-returns none but generates a conditional response matrix file. Let nS be the number of signals and nR be the
		number of responses the conditional response matrix is an nS X nR array with element i,j specifying the
		probability of the jth response to the ith input"""
	if lock_parameter is None:
		lock_parameter = {}
	with open(file_path, "r") as file:
		file_contents = file.readlines()
	header_string = ""
	
	if file_contents[1][0] != "#":
		column_index_list = ast.literal_eval(file_contents[0][2:])
	else:
		header_string += file_contents[0][2:-1]
		column_index_list = ast.literal_eval(file_contents[1][2:])
	moment_array = np.loadtxt(file_path, delimiter=",")
	# Formats single cell moment data into a 2d array

	if len(np.shape(moment_array)) == 1:
		moment_array = np.reshape(moment_array,(1,-1))
	if lock_parameter != {}:
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
	else:
		column_indices_to_keep = [i for i in range(np.shape(moment_array)[1])]
	moment_array = moment_array[:,column_indices_to_keep]
	crm = moment_array_to_crm(moment_array,assumed_distribution=assumed_distribution,discretization_parameter=discretization_parameter,percentile_cutoff=percentile_cutoff)
	np.savetxt(output_file_path,crm,header=header_string,delimiter=",")

