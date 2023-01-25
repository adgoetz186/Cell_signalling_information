import os
import numpy as np
from scipy.optimize import Bounds
from scipy.optimize import LinearConstraint
from scipy.optimize import minimize
import Mutual_Information_Final_Version.functions.data_processing_functions.Process_Raw_Cell_Conditioning as prc

def mutual_information_from_matrix(relative_input_vector, conditional_response_matrix):
    """obtains the mutual information from a conditional response matrix and a vector of input probabilities
    
    Input:
    
        -relative_input_vector (ndarray): Input distribution for signaling system
        
        -conditional_response_matrix (ndarray): Let nS be the number of signals and nR be the number of
        responses the conditional_probability_matrix is an nS X nR array with element i,j specifying the
        probability of the jth response to the ith input
        
    Output:
    
        -information_value: the mutual information of the specified conditional response matrix"""
    signal_ratios = [relative_input_vector]
    response_values = np.matmul(signal_ratios, conditional_response_matrix)
    signal_matrix = np.transpose(np.repeat(signal_ratios, repeats=np.shape(conditional_response_matrix)[1], axis=0))
    response_matrix = np.repeat(response_values, repeats=np.shape(conditional_response_matrix)[0], axis=0)
    ratio = np.divide(conditional_response_matrix, response_matrix, out=np.zeros_like(conditional_response_matrix), where=response_matrix != 0)
    log_of_ratio = np.log2(ratio, out=np.zeros_like(ratio), where=ratio != 0)
    mutual_info_matrix = np.multiply(np.multiply(signal_matrix, conditional_response_matrix), log_of_ratio)
    return -1*np.sum(np.nan_to_num(mutual_info_matrix, nan=0))


def mutual_information_from_crm_list_pp(relative_input_vector, list_of_crms, crm_index):
    """For use in parallel processing, finds the mutual information of one conditional response matrix from a list
    of conditional response matrices in a way that plays nice with parallel processing code
    
    Input:
    
        -relative_input_vector (ndarray): Input distribution for signaling system
        
        -crm_list (list): List of conditional response matrices
        
        -crm_index (int): specifies the index of the conditional response matrix
        within the "list_of_crms" for which the function will obtain the mutual information.
        
    Output:
    
        information_value: the mutual information of the specified conditional response matrix"""
    conditional_probability_matrix = list_of_crms[crm_index]
    return mutual_information_from_matrix(relative_input_vector,conditional_probability_matrix)


def conditional_mutual_information_from_list_of_crms(relative_input_vector, crm_list, return_average = True):
    """Calculates conditional mutual information from a list of CRMs
    
    Input:
    
        -relative_input_vector (ndarray): Input distribution for signaling system
        
        -crm_list (list): List of conditional response matrices
        
        -return_average (boolean): if true returns average information value of all CRMs,
        else gives list of information values
        
    Output:
    
        -information_value: list of each CRM information values or average of each CRM information values"""
    info_array = np.zeros(len(crm_list))
    for i in range(len(crm_list)):
        info_array[i] = mutual_information_from_matrix(relative_input_vector,crm_list[i])
    if return_average:
        return np.average(info_array)
    else:
        return info_array

def conditional_mutual_information_from_crm_folder(dir_or_file_path, input_distribution=None, return_all_cell_values=False, sort_filename=-1):
    """Returns the partially conditioned mutual information for the specified distribution of inputs

        Input:
        
            -dir_or_file_path (str): path to folder containing crm files or single crm file path
            
            -input_distribution (ndarray): distribution of inputs

            -return_all_cell_values (boolean): returns average of CRM performances or list of CRM performances

            -sort_filename (int): if -1 does not sort file, if positive gives the position of the sortable location in
            the file for example, a list of filenames with the ith file named "cell_number_i" would
            require "sort_filename = 2"

        Output:

            -information_value: either the average of the CRM mutual information values or a list of CRM mutual
            information values

            -input_dist: the channel capacity input distribution"""
    # Gets sorted list of conditional response matrices

    list_of_crm = []
    # Determines if path is directory or file
    if os.path.isdir(dir_or_file_path):
        file_name_list = os.listdir(dir_or_file_path)
        if sort_filename > -1:
            file_name_list = [file_name_list[i] for i in
                              np.argsort(np.array([int(i.split("_")[sort_filename]) for i in file_name_list]))]
        for crm_file in file_name_list:
            list_of_crm.append(np.loadtxt(dir_or_file_path + "/" + crm_file, delimiter=","))
    else:
        list_of_crm.append(np.loadtxt(dir_or_file_path, delimiter=","))
        
    # Assumes uniform input distribution if none specified
    if input_distribution is None:
        signal_count = np.shape(list_of_crm[0])[0]
        input_distribution = np.ones(signal_count) / signal_count

    # Obtains distribution of single cell performances under the given input distribution
    res = -1*conditional_mutual_information_from_list_of_crms(input_distribution, list_of_crm, not return_all_cell_values)
    return res
    
def pcmi_at_cc(dir_or_file_path, return_all_cell_values=False, sort_filename=-1):
    """Returns the partially conditioned mutual information at the channel capacity along with the channel capacity
    input distribution

    Input:
        -dir_or_file_path (str): path to folder containing crm files or single crm file path

        -return_all_cell_values (boolean): returns average of CRM performances or list of CRM performances

        -sort_filename (int): if -1 does not sort file, if positive gives the position of the sortable location in
        the file for example, a list of filenames with the ith file named "cell_number_i" would
        require "sort_filename = 2"

    Output:

        -information_value: either the average of the channel capacity of CRM mutual information values or a list of
        CRM mutual information values at the channel capacity

        -input_dist: the channel capacity input distribution"""
    # Gets sorted list of conditional response matrices
    list_of_crm = []
    # Determines if path is directory or file
    if os.path.isdir(dir_or_file_path):
        file_name_list = os.listdir(dir_or_file_path)
        if sort_filename > -1:
            file_name_list = [file_name_list[i] for i in
                              np.argsort(np.array([int(i.split("_")[sort_filename]) for i in file_name_list]))]
        for crm_file in file_name_list:
            list_of_crm.append(np.loadtxt(dir_or_file_path + "/" + crm_file, delimiter=","))
    else:
        list_of_crm.append(np.loadtxt(dir_or_file_path, delimiter=","))
    # Obtains channel capacity of list of conditional response matrices
    signal_count = np.shape(list_of_crm[0])[0]
    x0 = np.ones(signal_count) / signal_count
    bounds = Bounds(np.zeros(signal_count), np.ones(signal_count))
    linear_constraint = LinearConstraint(np.ones(signal_count), [1], [1])
    res = minimize(conditional_mutual_information_from_list_of_crms, x0, method='trust-constr',
                   constraints=linear_constraint, options={'verbose': 2}, bounds=bounds, args=(list_of_crm, True))
    # sns.heatmap(list_of_crm[0])
    # plt.show()
    # Returns results
    if return_all_cell_values:
        cell_mi = conditional_mutual_information_from_list_of_crms(res["x"], list_of_crm, False)
        return -1 * cell_mi, res["x"]
    else:
        return -1 * res["fun"], res["x"]


def pcmi_at_cc_generate_crm_on_fly(dir_path, lock_parameter=None, assumed_distribution="gamma",
                                   discretization_parameter=5, percentile_cutoff=0.0005):
    """Returns the partially conditioned mutual information at the channel capacity along with the channel capacity
	input distribution. Generates conditional response matrices without saving them.

	Input:

		-dir_path (str): path to folder containing moment files

		-lock_parameter (dict): specifies a dictionary of keys and values which must be in a column's dictionary
		for the moment to be used, often this will be a single entry dictionary of the form "{time:i_min}"

		-assumed_distribution (string): Specifies the assumed shape of the distribution of the conditional responses

		-discretization_parameter (int): Used to discretize continuous distributions using uniform binning,
		bin size is the inter quartile range of the smallest conditional response divided by the
		"discretization_parameter". It is approximately how many bins will be given to the IQR of the
		narrowest conditional response. Uniform binning is essential as entries within a column of the
		conditional response matrix MUST refer to the same responses.

		-percentile_cutoff (float): used only for discrete distributions to specify the effective min/max response
		by specifying the lowest percentile which should be included. Also removes the percentiles
		above "1-percentile_cutoff"

	Output:

		-information_value: either the average of the CRM mutual information values or a list of CRM mutual
		information values

		-input_dist: the channel capacity input distribution"""
    if lock_parameter is None:
        lock_parameter = {}
    # removes case dependence for easier comparison
    lc_assumed_dist = assumed_distribution.lower()
    # determines the number of moments required to describe the distribution and if the distribution is continuous
    continuous_distribution_moment_dictionary = {"gamma": 2, "gaussian": 2}
    discrete_distribution_moment_dictionary = {"poisson": 1, "nbinom": 2}
    if lc_assumed_dist in continuous_distribution_moment_dictionary.keys():
        is_cont = True
    else:
        is_cont = False
    if is_cont:
        number_of_moments = continuous_distribution_moment_dictionary[lc_assumed_dist]
    else:
        number_of_moments = discrete_distribution_moment_dictionary[lc_assumed_dist]
    
    # builds list of moments
    list_of_moments = []
    for moment_file in os.listdir(dir_path):
        moment_array, header = prc.moment_file_to_moments(dir_path, moment_file, lock_parameter=lock_parameter)
        list_of_moments.append(moment_array)
    # finds channel capacity of list of moments
    signal_count = int(np.shape(list_of_moments[0])[1] / number_of_moments)
    x0 = np.ones(signal_count) / signal_count
    bounds = Bounds(np.zeros(signal_count), np.ones(signal_count))
    linear_constraint = LinearConstraint(np.ones(signal_count), [1], [1])
    res = minimize(conditional_mutual_information_from_list_of_moments, x0, method='trust-constr',
                   constraints=linear_constraint, options={'verbose': 2}, bounds=bounds,
                   args=(list_of_moments, assumed_distribution, discretization_parameter, percentile_cutoff))
    return -1 * res["fun"], res["x"]

def conditional_mutual_information_from_list_of_moments(relative_input_vector, moment_array_list, assumed_distribution="gamma", discretization_parameter=5, percentile_cutoff=0.0005):
    """Calculates conditional mutual information from a list of moments
    
    Input:
        
        -relative_input_vector (ndarray): Input distribution for signaling system
        
        -moment_array_list (list): List of moment arrays, should be ordered with increments in dosage before
        increments in moment no. if applicable, (Moment 1 Dose 1, ..., Moment 1 Dose n, ..., Moment 2, Dose n)
        
        -assumed_distribution (string): Specifies the assumed shape of the distribution of the conditional responses
        
        -discretization_parameter (int): Used to discretize continuous distributions using uniform binning,
        bin size is the inter quartile range of the smallest conditional response divided by the
        "discretization_parameter". It is approximately how many bins will be given to the IQR of the
        narrowest conditional response. Uniform binning is essential as entries within a column of the
        conditional response matrix MUST refer to the same responses.
        
        -percentile_cutoff (float): used only for discrete distributions to specify the effective min/max response
        by specifying the lowest percentile which should be included. Also removes the percentiles
        above "1-percentile_cutoff"
        
    Output:
    
        -information_value: list of each CRM information values or average of each CRM information values"""
    info_array = np.zeros(len(moment_array_list))
    for i in range(len(moment_array_list)):
        # generates list of conditional response matrices
        conditional_probability_matrix = prc.moment_array_to_crm(moment_array_list[i],assumed_distribution = assumed_distribution,discretization_parameter=discretization_parameter, percentile_cutoff = percentile_cutoff)
        # uses list of conditional response matrices to get list of mutual information values
        info_array[i] = mutual_information_from_matrix(relative_input_vector,conditional_probability_matrix)
    return np.average(info_array)