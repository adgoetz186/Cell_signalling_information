import os
import numpy as np
import scipy.integrate as si
import mpmath as mp
from functools import partial
import scipy.stats as st
from scipy.optimize import Bounds
from scipy.optimize import LinearConstraint
from scipy.optimize import minimize
import Mutual_Information_Main.functions.data_processing_functions.Process_Raw_Cell_Conditioning as prc

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


def pcmi_at_cc_from_crm_list(list_of_crms, return_all_cell_values=False):
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
    
    # finds channel capacity of list of moments
    signal_count = int(np.shape(list_of_crms[0])[0])
    x0 = np.ones(signal_count) / signal_count
    bounds = Bounds(np.zeros(signal_count), np.ones(signal_count))
    linear_constraint = LinearConstraint(np.ones(signal_count), [1], [1])
    res = minimize(conditional_mutual_information_from_list_of_crms, x0, method='trust-constr',
                   constraints=linear_constraint, options={'verbose': 2}, bounds=bounds,
                   args=(list_of_crms))
    if return_all_cell_values:
        cell_mi = conditional_mutual_information_from_list_of_crms(res["x"], list_of_crms, return_average = False )
        return -1 * cell_mi, res["x"]
    else:
        return -1 * res["fun"], res["x"]

def integrand_for_direct_integration(r,relative_input_vector,moment_array):
    prgu_array = np.zeros_like(relative_input_vector)
    means = moment_array[0,:np.size(relative_input_vector)]
    vars = moment_array[0,np.size(relative_input_vector):] - means**2
    scales = vars / means
    shapes = means / scales
    for u_ind in range(np.size(relative_input_vector)):
        prgu_array[u_ind] = st.gamma.pdf(r,shapes[u_ind],scale=scales[u_ind])
    denom = np.sum(prgu_array*relative_input_vector)
    return np.sum(np.nan_to_num(prgu_array*relative_input_vector*np.log2(prgu_array/denom),nan=0))
        
    

def mutual_information_direct_integration_gamma(relative_input_vector, moment_array, percentile_cutoff=0.000005,error=1e-7):
    """Calculates conditional mutual information from a list of moment array. Numerically integrates gamma distribution
        rather than binning

    Input:

        -relative_input_vector (ndarray): Input distribution for signaling system

        -moment_array (ndarray): moment array, should be ordered with increments in dosage before
        increments in moment no. if applicable, (Moment 1 Dose 1, ..., Moment 1 Dose n, ..., Moment 2, Dose n)
        
        -response_bounds (ndarray): min and max values for the response to use in integration

    Output:

        -information_value: information for the specified system"""
    means = moment_array[0, :np.size(relative_input_vector)]
    vars = moment_array[0, np.size(relative_input_vector):] - means ** 2
    scales = vars / means
    shapes = means / scales
    min_list = []
    max_list = []
    for u_ind in range(np.size(relative_input_vector)):
        min_list.append(st.gamma.ppf(percentile_cutoff, shapes[u_ind], scale=scales[u_ind]))
        max_list.append(st.gamma.ppf(1 - percentile_cutoff, shapes[u_ind], scale=scales[u_ind]))
    return si.quad(integrand_for_direct_integration,np.min(np.array(min_list)),np.max(np.array(max_list)),args=(relative_input_vector,moment_array),limit=1000,epsabs=error,epsrel = error)