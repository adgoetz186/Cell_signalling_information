import numpy as np
import scipy.stats as st
import scipy.special as sp
from Mutual_Information_Final_Version.functions.mutual_information_functions.MI_Calculation_Functions import mutual_information_from_matrix


def pxgut(ulist,kdeg,rec_val,kbind,kunbind,max_response_factor):
    """Generates a conditional response matrix for a given parameter set.
        The matrices' rows correspond to inputs (u) and the columns correspond to responses.

            Input:

                -ulist (ndarray): List of inputs

                -kdeg (flaot): Value of k_deg

                -rec_val (int): number of average receptors at t = 0

                -kbind (float): binding rate of ligand

                -kunbind (float): unbinding rate of ligand

                -max_response_factor (int): sets max number of ligand bound receptors allowed for the system based on mean
                max allowed ligand bound receptors = max_response_factor*rec_val

            Output:

                -list_of_crm (list): A list containing a single CRM for the specified parameters, a list is used
                to mirror other pxgut functions"""
    # Declares the list of allowed responses
    xlist = np.arange(0, rec_val*max_response_factor)

    # generates list of mean response values for each cell and response
    mean_array = np.empty((np.size(ulist),1))
    for uind in range(np.size(ulist)):
        mean_array[uind] = ulist[uind] * rec_val * kbind / (ulist[uind] * kbind + kdeg + kunbind)

    # generates crm of the cell with the specified parameters
    pxgutdist = np.empty((np.size(ulist),np.size(xlist)))
    for uind in range(np.size(ulist)):
        # the distribution of a cell's response to a given input
        pxgutdist[uind] = st.poisson.pmf(xlist, mean_array[uind, 0])
        # Ensures the probability distribution sums to 1 to safeguard against numerical issues
        pxgutdist[uind] /= np.sum(pxgutdist[uind])
    return pxgutdist


def pxgut_deg(ulist, tvar, tmean, paradraws, rec_val, kbind, kunbind, max_response_factor):
    """Generates a list of conditional response matrices, each CRM represents a unique cell parameter value (t).
    The matrices' rows correspond to inputs (u) and the columns correspond to responses.
    Assumes a cell population in which kdeg is the only varying parameter.

        Input:

            -ulist (ndarray): List of inputs

            -tvar (flaot): Variance of k_deg
            
            -tmean (flaot): Mean of k_deg
            
            -paradraws (int): number of parameters to draw (or the number of cells to simulate)
            
            -rec_val (int): number of average receptors at t = 0
            
            -kbind (float): binding rate of ligand
            
            -kunbind (float): unbinding rate of ligand
            
            -max_response_factor (int): sets max number of ligand bound receptors allowed for the system based on mean
            max allowed ligand bound receptors = max_response_factor*rec_val

        Output:

            -list_of_crm (list): A list of conditional response matrices, each matrix belonging to a simulated cell
            with its own value of kdeg"""
    
    # Declares the list of allowed responses
    xlist = np.arange(0, rec_val*max_response_factor)
    
    # generates list of theta values
    theta_scale1 = tvar / tmean
    theta_shape1 = tmean / theta_scale1
    kdeglist = np.random.gamma(shape=theta_shape1, scale=theta_scale1, size=(paradraws))
    
    # generates list of mean response values for each cell and response
    mean_array = np.empty((np.size(ulist),np.size(kdeglist)))
    for uind in range(np.size(ulist)):
        mean_array[uind] = ulist[uind] * rec_val * kbind / (ulist[uind] * kbind + kdeglist + kunbind)
        
    # generates list of crms where each crm represents the crm of a randomly generated cell
    list_of_crm = []
    for tind in range(len(kdeglist)):
        
        pxgutdist = np.empty((np.size(ulist),np.size(xlist)))
        for uind in range(np.size(ulist)):
            # the distribution of a cell's response to a given input
            pxgutdist[uind] = st.poisson.pmf(xlist, mean_array[uind, tind])
            # Ensures the probability distribution sums to 1 to safeguard against numerical issues
            pxgutdist[uind] /= np.sum(pxgutdist[uind])
        list_of_crm.append(pxgutdist)
    return list_of_crm


def pxgut_r0(ulist,tvar,tmean,paradraws,deg_val,kbind,kunbind,max_response_factor):
    """Generates a list of conditional response matrices, each CRM represents a unique cell parameter value (t).
        The matrices' rows correspond to inputs (u) and the columns correspond to responses.
        Assumes a cell population in which r0 is the only varying parameter.
        
            Input:

                -ulist (ndarray): List of inputs

                -tvar (flaot): Variance of r_0

                -tmean (flaot): Mean of r_0

                -paradraws (int): number of parameters to draw (or the number of cells to simulate)

                -deg_val (float): value of k_deg

                -kbind (float): binding rate of ligand

                -kunbind (float): unbinding rate of ligand

                -max_response_factor (int): sets max number of ligand bound receptors allowed for the system based on
                mean; max allowed ligand bound receptors = max_response_factor*tmean

            Output:

                -list_of_crm (list): A list of conditional response matrices, each matrix belonging to a simulated cell
                with its own value of kdeg"""
    # Declares the list of allowed responses
    xlist = np.arange(0, tmean*max_response_factor)

    # generates list of theta values
    theta_scale1 = tvar / tmean
    theta_shape1 = tmean / theta_scale1
    r0list = np.random.gamma(shape=theta_shape1, scale=theta_scale1, size=(paradraws))

    # generates list of mean response values for each cell and response
    mean_array = np.empty((np.size(ulist),np.size(r0list)))
    for uind in range(np.size(ulist)):
        mean_array[uind] = ulist[uind] * r0list * kbind / (ulist[uind] * kbind + deg_val + kunbind)

    # generates list of crms where each crm represents the crm of a randomly generated cell
    list_of_crm = []
    for tind in range(len(r0list)):
        pxgutdist = np.empty((np.size(ulist),np.size(xlist)))
        for uind in range(np.size(ulist)):
            # the distribution of a cell's response to a given input
            pxgutdist[uind] = st.poisson.pmf(xlist,mean_array[uind,tind])
            
            # Ensures the probability distribution sums to 1 to safeguard against numerical issues
            pxgutdist[uind] /= np.sum(pxgutdist[uind])
            
        list_of_crm.append(pxgutdist)
    return list_of_crm


def single_cell_MI(ulist,kdeg,rec_val,kbind,kunbind,max_response_factor = 5):
    """Generates a single cell state's mutual information value assumes a uniform distribution over inputs.

            Input:

                -ulist (ndarray): List of inputs

                -kdeg (flaot): Value of k_deg

                -rec_val (int): number of average receptors at t = 0

                -kbind (float): binding rate of ligand

                -kunbind (float): unbinding rate of ligand

                -max_response_factor (int): sets max number of ligand bound receptors allowed for the system based on mean
                max allowed ligand bound receptors = max_response_factor*rec_val

            Output:

                -single_cell_MI (float): The single cell state's information value"""
    # generates list of conditional response matrices for a randomly generated collection of cell states
    crm = pxgut(ulist,kdeg,rec_val,kbind,kunbind,max_response_factor)

    # generates input signal array
    signal_dis = np.ones(np.size(ulist))
    signal_dis /= np.sum(signal_dis)

    # returns the mutual information for the cell
    return mutual_information_from_matrix(signal_dis,crm)


def cmi_deg_var(ulist,tvar,tmean,paradraws,rec_val,kbind,kunbind,max_response_factor = 5):
    """Generates a list of single cell state's mutual information values assumes a uniform distribution over inputs.
    Assumes a cell population in which kdeg is the only varying parameter.

            Input:

                -ulist (ndarray): List of inputs

                -tvar (flaot): Variance of k_deg

                -tmean (flaot): Mean of k_deg

                -paradraws (int): number of parameters to draw (or the number of cells/cell states to simulate)

                -rec_val (int): number of average receptors at t = 0

                -kbind (float): binding rate of ligand

                -kunbind (float): unbinding rate of ligand

                -max_response_factor (int): sets max number of ligand bound receptors allowed for the system based on mean
                max allowed ligand bound receptors = max_response_factor*rec_val

            Output:

                -single_cell_MI (ndarray): An array with each entry being a single cell state's information value"""
    # generates list of conditional response matrices for a randomly generated collection of cell states
    crm_list = pxgut_deg(ulist,tvar,tmean,paradraws,rec_val,kbind,kunbind,max_response_factor)
    
    # generates input signal array
    signal_dis = np.ones(np.size(ulist))
    signal_dis /= np.sum(signal_dis)
    
    # obtains the mutual information for each cell
    single_cell_MI = np.zeros((len(crm_list)))
    for crm_ind in range(len(crm_list)):
        single_cell_MI[crm_ind] = mutual_information_from_matrix(signal_dis, crm_list[crm_ind])
    return single_cell_MI


def cmi_r0_var(ulist,tvar,tmean,paradraws,deg_val,kbind,kunbind,max_response_factor = 5):
    """Generates a list of single cell state's mutual information values assumes a uniform distribution over inputs.
        Assumes a cell population in which r0 is the only varying parameter.

            Input:

                -ulist (ndarray): List of inputs

                -tvar (flaot): Variance of r_0

                -tmean (flaot): Mean of r_0

                -paradraws (int): number of parameters to draw (or the number of cells/cell states to simulate)

                -deg_val (int): value of k_deg

                -kbind (float): binding rate of ligand

                -kunbind (float): unbinding rate of ligand

                -max_response_factor (int): sets max number of ligand bound receptors allowed for the system based
                on mean; max allowed ligand bound receptors = max_response_factor*tmean

            Output:

                -single_cell_MI (ndarray): An array with each entry being a single cell state's information value"""
    # generates list of conditional response matrices for a randomly generated collection of cell states
    crm_list = pxgut_r0(ulist,tvar,tmean,paradraws,deg_val,kbind,kunbind,max_response_factor)

    # generates input signal array
    signal_dis = np.ones(np.size(ulist))
    signal_dis /= np.sum(signal_dis)

    # obtains the mutual information for each cell
    single_cell_MI = np.zeros((len(crm_list)))
    for crm_ind in range(len(crm_list)):
        single_cell_MI[crm_ind] = mutual_information_from_matrix(signal_dis, crm_list[crm_ind])
    return single_cell_MI


def pxgu_deg(u,tvar,tmean,paradraws,rec_val,kbind,kunbind,max_response_factor):
    """Generates a list of single cell state's mutual information values assumes a uniform distribution over inputs.
        Assumes a cell population in which k_deg is the only varying parameter.

            Input:

                -u (float): input value

                -tvar (flaot): Variance of k_deg

                -tmean (flaot): Mean of k_deg

                -paradraws (int): number of parameters to draw (or the number of cells/cell states to simulate)

                -rec_val (int): value of r_0

                -kbind (float): binding rate of ligand

                -kunbind (float): unbinding rate of ligand
                
                -max_response_factor (int): sets max number of ligand bound receptors allowed for the system based
                on mean; max allowed ligand bound receptors = max_response_factor*rec_val

            Output:

                -pxgutdist (ndarray): The conditional response to a given input for a given cell state"""
    # Declares the list of allowed responses
    xlist = np.arange(0, rec_val * max_response_factor)
    
    # generates list of theta values
    theta_scale1 = tvar / tmean
    theta_shape1 = tmean / theta_scale1
    kdeglist = np.random.gamma(shape=theta_shape1, scale=theta_scale1, size=(paradraws))

    # generates list of mean response values for each cell and response
    mean_array = u*rec_val*kbind/(u*kbind+kdeglist+kunbind)
    
    # finds the first moment of the population
    one_mom_pmf = np.average(mean_array)
    # finds the second moment of the population
    two_mom_pmf = np.average(mean_array**2+mean_array)
    # finds the variance of the population
    var_pmf = two_mom_pmf - one_mom_pmf**2
    
    # obtains the distribution of the population assuming a negative binomial distribution
    p = one_mom_pmf/var_pmf
    n = one_mom_pmf**2/(var_pmf-one_mom_pmf)
    pxgutdist = st.nbinom.pmf(xlist,n,p)
    
    # Ensures the probability distribution sums to 1 to safeguard against numerical issues
    pxgutdist /= np.sum(pxgutdist)
    return pxgutdist


def pxgu_rec(u,tvar,tmean,paradraws,degvalue,kbind,kunbind,max_response_factor):
    """Generates a list of single cell state's mutual information values assumes a uniform distribution over inputs.
        Assumes a cell population in which r_0 is the only varying parameter.

            Input:

                -u (float): input value

                -tvar (flaot): Variance of r_0

                -tmean (flaot): Mean of r_0

                -paradraws (int): number of parameters to draw (or the number of cells/cell states to simulate)

                -degvalue (int): value of k_deg

                -kbind (float): binding rate of ligand

                -kunbind (float): unbinding rate of ligand

                -max_response_factor (int): sets max number of ligand bound receptors allowed for the system based
                on mean; max allowed ligand bound receptors = max_response_factor*tmean

            Output:

                -pxgutdist (ndarray): The conditional response to a given input for a given cell state"""
    # Declares the list of allowed responses
    xlist = np.arange(0, tmean * max_response_factor)

    # generates list of theta values
    theta_scale1 = tvar / tmean
    theta_shape1 = tmean / theta_scale1
    R0list = np.random.gamma(shape=theta_shape1, scale=theta_scale1, size=(paradraws))

    # generates list of mean response values for each cell and response
    mean_array = u*R0list*kbind/(u*kbind+degvalue+kunbind)

    # finds the first moment of the population
    one_mom_pmf = np.average(mean_array)

    # finds the second moment of the population
    two_mom_pmf = np.average(mean_array ** 2 + mean_array)

    # finds the variance of the population
    var_pmf = two_mom_pmf - one_mom_pmf ** 2

    # obtains the distribution of the population assuming a negative binomial distribution
    p = one_mom_pmf / var_pmf
    n = one_mom_pmf ** 2 / (var_pmf - one_mom_pmf)
    pxgutdist = st.nbinom.pmf(xlist, n, p)

    # Ensures the probability distribution sums to 1 to safeguard against numerical issues
    pxgutdist /= np.sum(pxgutdist)
    return pxgutdist


def mi_deg_var(ulist,tvar,tmean,paradraws,rec_val,kbind,kunbind,max_response_factor = 5):
    """Generates the MI of the CSAR for the toy model population, assumes a uniform distribution over inputs.
        Assumes a cell population in which k_deg is the only varying parameter.

            Input:

                -ulist (ndarray): List of inputs

                -tvar (flaot): Variance of k_deg

                -tmean (flaot): Mean of k_deg

                -paradraws (int): number of parameters to draw (or the number of cells/cell states to simulate)

                -rec_val (int): value of r_0

                -kbind (float): binding rate of ligand

                -kunbind (float): unbinding rate of ligand

                -max_response_factor (int): sets max number of ligand bound receptors allowed for the system based
                on mean; max allowed ligand bound receptors = max_response_factor*tmean

            Output:

                -mi_of_csar (float): The mutual information of the cell state average response"""
    # generates conditional response matrix for cell state averaged responses
    crm = np.zeros((np.size(ulist),max_response_factor*rec_val))
    for uind in range(np.size(ulist)):
        crm[uind] = pxgu_deg(ulist[uind],tvar,tmean,paradraws,rec_val,kbind,kunbind,max_response_factor)
        
    # creates uniform signal distribution
    signal_dis = np.ones(np.size(ulist))/np.size(ulist)
    
    # obtains the mutual information
    mi_of_csar = mutual_information_from_matrix(signal_dis, crm)
    return mi_of_csar


def mi_r0_var(ulist,tvar,tmean,paradraws,deg_val,kbind,kunbind,max_response_factor = 5):
    """Generates the MI of the CSAR for the toy model population, assumes a uniform distribution over inputs.
        Assumes a cell population in which r_0 is the only varying parameter.

            Input:

                -ulist (ndarray): List of inputs

                -tvar (flaot): Variance of r_0

                -tmean (flaot): Mean of r_0

                -paradraws (int): number of parameters to draw (or the number of cells/cell states to simulate)

                -deg_val (int): value of k_deg

                -kbind (float): binding rate of ligand

                -kunbind (float): unbinding rate of ligand

                -max_response_factor (int): sets max number of ligand bound receptors allowed for the system based
                on mean; max allowed ligand bound receptors = max_response_factor*tmean

            Output:

                -mi_of_csar (float): The mutual information of the cell state average response"""
    # generates conditional response matrix for cell state averaged responses
    crm = np.zeros((np.size(ulist),tmean * max_response_factor))
    for uind in range(np.size(ulist)):
        crm[uind] = pxgu_rec(ulist[uind],tvar,tmean,paradraws,deg_val,kbind,kunbind,max_response_factor)

    # creates uniform signal distribution
    signal_dis = np.ones(np.size(ulist))/np.size(ulist)

    # obtains the mutual information
    mi_of_csar = mutual_information_from_matrix(signal_dis, crm)
    return mi_of_csar


def build_heatmap_from_sc_database(sc_db):
    """Takes the single cell pandas database used to store outputs and generates a ndarray which can be easily
    plotted.

        Input:

            -sc_db (database): pandas database which stores population performance values

        Output:

            -heatmap_array (ndarray): Matrix of information values where element aij represents the mutual information
             for the ith input distribution and the jth parameter distribution."""
    input_list = list(set(sc_db["input_cv"]))
    input_list.sort()
    param_list = list(set(sc_db["parameter_cv"]))
    param_list.sort()
    heatmap_array = np.zeros((np.size(input_list),np.size(param_list)))
    for input_ind in range(len(input_list)):
        for param_ind in range(len(param_list)):
            df_ss = sc_db[(sc_db["input_cv"] == input_list[input_ind]) & (sc_db["parameter_cv"] == param_list[param_ind])]
            heatmap_array[input_ind,param_ind] = np.average(df_ss["cell_mi"].to_numpy())
    return heatmap_array
            