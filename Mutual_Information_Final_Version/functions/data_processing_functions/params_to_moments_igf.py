import numpy as np
from scipy.stats import gamma
from scipy.integrate import solve_ivp

# This code was written by Hoda Akl with moment differential equations obtained by Andrew Goetz

def Get_Init_fn(k):
    """Inputs: k
    Output: z0 """
    # last two entries
    nCom = 44
    k = 10**(k)

    k1 = k[0]*k[1]   #Synthesis of IGFR
    k2 = k[1]   #Degredation of IGFR
    # k3 = k[2]   #Binding of IGFR to IGF
    # k4 = k[3]   #Unbinding IGFR to IGF
    # k5 = k[4]   #Phosphorylation of bound receptor
    # k6 = k[5]   #Dephosphorylation of bound receptor
    # k7 = k[6]   #Dephosphorylation of AKT
    # k8 = k[7]   #Phosphorylation of AKT
    # k9 = k[8]   #Phosphorylation of FoxO
    # k10 = k[9]  #Dephosphorylation of FoxO
    k11 = k[10] #Influx of FoxO to nucleus
    k12 = k[11] #Efflux of FoxO from nucleus
    k_tot_Akt = k[12]
    k_tot_foxo = k[13]

    #
    z0 = np.zeros(nCom)
    # z0_1
    z0[0]= k1/k2
    # z0_3
    z0[3] = k_tot_Akt # 34050   # total AKT which can vary
    # z0_7
    z0[6] = (k_tot_foxo*k12)/(k11 + k12)   #710 is total FoxO which can also vary
    # z0_8
    z0[7] = (k_tot_foxo*k11)/(k11 + k12)
    # z0_1_1
    z0[8]= k1/k2
    # z0_7_7
    # -(- 710*k11^2*k12 + 710*k11*k12^2)/(k11^2*(k11 + k12))
    z0[41] = (k_tot_foxo*k11*k12)/(k11**2 + 2*k11*k12 + k12**2)
    # z0_7_8
    z0[42] =  -(k_tot_foxo*k11*k12)/((k11 + k12)**2)
    # z0_8_8
    z0[43] = (k_tot_foxo*k11*k12)/(k11**2 + 2*k11*k12 + k12**2)
    return z0

def MomentsDiff_Eq_fn(t, z, k, IGF):
	"""Inputs: t (in seconds), z(array of 44) , k (array of len 12), IGF (concentration in pM)
	Outputs: a list of differential equations """
	# make the exponent of the log rates
	# define the rates
	
	k = 10 ** (k)
	
	k1 = k[0] * k[1]  # number of receptors
	k2 = k[1]  # Degredation of IGFR
	k3 = k[2]  # Binding of IGFR to IGF
	k4 = k[3]  # Unbinding IGFR to IGF
	k5 = k[4]  # Phosphorylation of bound receptor
	k6 = k[5]  # Dephosphorylation of bound receptor
	k7 = k[6]  # Dephosphorylation of AKT
	k8 = k[7]  # Phosphorylation of AKT
	k9 = k[8]  # Phosphorylation of FoxO
	k10 = k[9]  # Dephosphorylation of FoxO
	k11 = k[10]  # Influx of FoxO to nucleus
	k12 = k[11]  # Efflux of FoxO from nucleus
	# define the species
	# z_1 = z[0]   #R
	# z_2 = z[1]   #B
	# z_3 = z[2]   #P
	# z_4 = z[3]   #akt
	# z5 = z[4]   #pakt
	# z6 = z[5]   #pfoxoc
	# z7 = z[6]   #foxoc
	# z8 = z[7]   #foxon
	# zi is the mean of i
	# zii is the 2nd moment of i
	# zij is the comoment of i & j
	z_1, z_2, z_3, z_4, z_5, z_6, z_7, z_8 = z[:8]
	z_1_1, z_1_2, z_1_3, z_1_4, z_1_5, z_1_6, z_1_7, z_1_8 = z[8:16]
	z_2_2, z_2_3, z_2_4, z_2_5, z_2_6, z_2_7, z_2_8 = z[16:23]
	z_3_3, z_3_4, z_3_5, z_3_6, z_3_7, z_3_8 = z[23:29]
	z_4_4, z_4_5, z_4_6, z_4_7, z_4_8 = z[29:34]
	z_5_5, z_5_6, z_5_7, z_5_8 = z[34:38]
	z_6_6, z_6_7, z_6_8 = z[38:41]
	z_7_7, z_7_8 = z[41:43]
	z_8_8 = z[43]
	
	# initialize the differential equation array
	ns = 44
	dsyn = np.zeros(ns)
	# write the differential equations
	# z_1'
	dsyn[0] = k1 - k2 * z_1 - IGF * k3 * z_1 + k4 * z_2
	# z_2'
	dsyn[1] = IGF * k3 * z_1 - k2 * z_2 - k4 * z_2 - k5 * z_2 + k6 * z_3
	# z_3'
	dsyn[2] = k5 * z_2 - k2 * z_3 - k6 * z_3
	# z_4'
	dsyn[3] = -k8 * z_3 * z_4 + k7 * z_5 - k8 * z_3_4
	# z_5'
	dsyn[4] = k8 * z_3 * z_4 - k7 * z_5 + k8 * z_3_4
	# z_6'
	dsyn[5] = -k10 * z_6 + k9 * z_5 * z_7 + k9 * z_5_7
	# z_7'
	dsyn[6] = k10 * z_6 - k11 * z_7 - k9 * z_5 * z_7 + k12 * z_8 - k9 * z_5_7
	# z_8'
	dsyn[7] = k11 * z_7 - k12 * z_8
	# z_1_1'
	dsyn[8] = k1 + k2 * z_1 + IGF * k3 * z_1 + k4 * z_2 - 2 * k2 * z_1_1 - 2 * IGF * k3 * z_1_1 + 2 * k4 * z_1_2
	# z_1_2'
	dsyn[
		9] = -IGF * k3 * z_1 - k4 * z_2 + IGF * k3 * z_1_1 - 2 * k2 * z_1_2 - IGF * k3 * z_1_2 - k4 * z_1_2 - k5 * z_1_2 + k6 * z_1_3 + k4 * z_2_2
	# z_1_3'
	dsyn[10] = k5 * z_1_2 - 2 * k2 * z_1_3 - IGF * k3 * z_1_3 - k6 * z_1_3 + k4 * z_2_3
	# z_1_4'
	dsyn[11] = -k8 * z_4 * z_1_3 - k2 * z_1_4 - IGF * k3 * z_1_4 - k8 * z_3 * z_1_4 + k7 * z_1_5 + k4 * z_2_4
	# z_1_5'
	dsyn[12] = k8 * z_4 * z_1_3 + k8 * z_3 * z_1_4 - k2 * z_1_5 - IGF * k3 * z_1_5 - k7 * z_1_5 + k4 * z_2_5
	# z_1_6'
	dsyn[13] = k9 * z_7 * z_1_5 - k10 * z_1_6 - k2 * z_1_6 - IGF * k3 * z_1_6 + k9 * z_5 * z_1_7 + k4 * z_2_6
	# z_1_7'
	dsyn[
		14] = -k9 * z_7 * z_1_5 + k10 * z_1_6 - k11 * z_1_7 - k2 * z_1_7 - IGF * k3 * z_1_7 - k9 * z_5 * z_1_7 + k12 * z_1_8 + k4 * z_2_7
	# z_1_8'
	dsyn[15] = k11 * z_1_7 - k12 * z_1_8 - k2 * z_1_8 - IGF * k3 * z_1_8 + k4 * z_2_8
	# z_2_2'
	dsyn[
		16] = IGF * k3 * z_1 + k2 * z_2 + k4 * z_2 + k5 * z_2 + k6 * z_3 + 2 * IGF * k3 * z_1_2 - 2 * k2 * z_2_2 - 2 * k4 * z_2_2 - 2 * k5 * z_2_2 + 2 * k6 * z_2_3
	# z_2_3'
	dsyn[
		17] = -k5 * z_2 - k6 * z_3 + IGF * k3 * z_1_3 + k5 * z_2_2 - 2 * k2 * z_2_3 - k4 * z_2_3 - k5 * z_2_3 - k6 * z_2_3 + k6 * z_3_3
	# z_2_4'
	dsyn[
		18] = IGF * k3 * z_1_4 - k8 * z_4 * z_2_3 - k2 * z_2_4 - k4 * z_2_4 - k5 * z_2_4 - k8 * z_3 * z_2_4 + k7 * z_2_5 + k6 * z_3_4
	# z_2_5'
	dsyn[
		19] = IGF * k3 * z_1_5 + k8 * z_4 * z_2_3 + k8 * z_3 * z_2_4 - k2 * z_2_5 - k4 * z_2_5 - k5 * z_2_5 - k7 * z_2_5 + k6 * z_3_5
	# z_2_6'
	dsyn[
		20] = IGF * k3 * z_1_6 + k9 * z_7 * z_2_5 - k10 * z_2_6 - k2 * z_2_6 - k4 * z_2_6 - k5 * z_2_6 + k9 * z_5 * z_2_7 + k6 * z_3_6
	# z_2_7'
	dsyn[
		21] = IGF * k3 * z_1_7 - k9 * z_7 * z_2_5 + k10 * z_2_6 - k11 * z_2_7 - k2 * z_2_7 - k4 * z_2_7 - k5 * z_2_7 - k9 * z_5 * z_2_7 + k12 * z_2_8 + k6 * z_3_7
	# z_2_8'
	dsyn[22] = IGF * k3 * z_1_8 + k11 * z_2_7 - k12 * z_2_8 - k2 * z_2_8 - k4 * z_2_8 - k5 * z_2_8 + k6 * z_3_8
	# z_3_3'
	dsyn[23] = k5 * z_2 + k2 * z_3 + k6 * z_3 + 2 * k5 * z_2_3 - 2 * k2 * z_3_3 - 2 * k6 * z_3_3
	# z_3_4'
	dsyn[24] = k5 * z_2_4 - k8 * z_4 * z_3_3 - k2 * z_3_4 - k6 * z_3_4 - k8 * z_3 * z_3_4 + k7 * z_3_5
	# z_3_5'
	dsyn[25] = k5 * z_2_5 + k8 * z_4 * z_3_3 + k8 * z_3 * z_3_4 - k2 * z_3_5 - k6 * z_3_5 - k7 * z_3_5
	# z_3_6'
	dsyn[26] = k5 * z_2_6 + k9 * z_7 * z_3_5 - k10 * z_3_6 - k2 * z_3_6 - k6 * z_3_6 + k9 * z_5 * z_3_7
	# z_3_7'
	dsyn[
		27] = k5 * z_2_7 - k9 * z_7 * z_3_5 + k10 * z_3_6 - k11 * z_3_7 - k2 * z_3_7 - k6 * z_3_7 - k9 * z_5 * z_3_7 + k12 * z_3_8
	# z_3_8'
	dsyn[28] = k5 * z_2_8 + k11 * z_3_7 - k12 * z_3_8 - k2 * z_3_8 - k6 * z_3_8
	# z_4_4'
	dsyn[29] = k8 * z_3 * z_4 + k7 * z_5 + k8 * z_3_4 - 2 * k8 * z_4 * z_3_4 - 2 * k8 * z_3 * z_4_4 + 2 * k7 * z_4_5
	# z_4_5'
	dsyn[
		30] = -k8 * z_3 * z_4 - k7 * z_5 - k8 * z_3_4 + k8 * z_4 * z_3_4 - k8 * z_4 * z_3_5 + k8 * z_3 * z_4_4 - k7 * z_4_5 - k8 * z_3 * z_4_5 + k7 * z_5_5
	# z_4_6'
	dsyn[31] = -k8 * z_4 * z_3_6 + k9 * z_7 * z_4_5 - k10 * z_4_6 - k8 * z_3 * z_4_6 + k9 * z_5 * z_4_7 + k7 * z_5_6
	# z_4_7'
	dsyn[
		32] = -k8 * z_4 * z_3_7 - k9 * z_7 * z_4_5 + k10 * z_4_6 - k11 * z_4_7 - k8 * z_3 * z_4_7 - k9 * z_5 * z_4_7 + k12 * z_4_8 + k7 * z_5_7
	# z_4_8'
	dsyn[33] = -k8 * z_4 * z_3_8 + k11 * z_4_7 - k12 * z_4_8 - k8 * z_3 * z_4_8 + k7 * z_5_8
	# z_5_5'
	dsyn[34] = k8 * z_3 * z_4 + k7 * z_5 + k8 * z_3_4 + 2 * k8 * z_4 * z_3_5 + 2 * k8 * z_3 * z_4_5 - 2 * k7 * z_5_5
	# z_5_6'
	dsyn[35] = k8 * z_4 * z_3_6 + k8 * z_3 * z_4_6 + k9 * z_7 * z_5_5 - k10 * z_5_6 - k7 * z_5_6 + k9 * z_5 * z_5_7
	# z_5_7'
	dsyn[
		36] = k8 * z_4 * z_3_7 + k8 * z_3 * z_4_7 - k9 * z_7 * z_5_5 + k10 * z_5_6 - k11 * z_5_7 - k7 * z_5_7 - k9 * z_5 * z_5_7 + k12 * z_5_8
	# z_5_8'
	dsyn[37] = k8 * z_4 * z_3_8 + k8 * z_3 * z_4_8 + k11 * z_5_7 - k12 * z_5_8 - k7 * z_5_8
	# z_6_6'
	dsyn[38] = k10 * z_6 + k9 * z_5 * z_7 + 2 * k9 * z_7 * z_5_6 + k9 * z_5_7 - 2 * k10 * z_6_6 + 2 * k9 * z_5 * z_6_7
	# z_6_7'
	dsyn[
		39] = -k10 * z_6 - k9 * z_5 * z_7 - k9 * z_7 * z_5_6 - k9 * z_5_7 + k9 * z_7 * z_5_7 + k10 * z_6_6 - k10 * z_6_7 - k11 * z_6_7 - k9 * z_5 * z_6_7 + k12 * z_6_8 + k9 * z_5 * z_7_7
	# z_6_8'
	dsyn[40] = k9 * z_7 * z_5_8 + k11 * z_6_7 - k10 * z_6_8 - k12 * z_6_8 + k9 * z_5 * z_7_8
	# z_7_7'
	dsyn[
		41] = k10 * z_6 + k11 * z_7 + k9 * z_5 * z_7 + k12 * z_8 + k9 * z_5_7 - 2 * k9 * z_7 * z_5_7 + 2 * k10 * z_6_7 - 2 * k11 * z_7_7 - 2 * k9 * z_5 * z_7_7 + 2 * k12 * z_7_8
	# z_7_8'
	dsyn[
		42] = -k11 * z_7 - k12 * z_8 - k9 * z_7 * z_5_8 + k10 * z_6_8 + k11 * z_7_7 - k11 * z_7_8 - k12 * z_7_8 - k9 * z_5 * z_7_8 + k12 * z_8_8
	# z_8_8'
	dsyn[43] = k11 * z_7 + k12 * z_8 + 2 * k11 * z_7_8 - 2 * k12 * z_8_8
	#
	return dsyn

def solve_Moments_fn(K,IGF,t, z0, meth = 'BDF'):
    """Inputs: K , L = IGF , tend
    Outputs: Z solution"""
    tspan = [0, t]
    # z0 = Get_Init_fn(K)
    sol_dyn = solve_ivp(MomentsDiff_Eq_fn, tspan, z0, method = meth, args=(K,IGF))
    sol = sol_dyn.y[:,-1]
    return sol_dyn, sol

def FoxOn_preds_fn(k, L, t, z0, meth='BDF'):
	""""Function spits out mean and second moment of nuclear foxO at time t
	input:  L = IGF concentration in nM
			t = end time in seconds
	output: meanfoxO, variancefoxO"""
	_, Sol_t = solve_Moments_fn(k, L, t, z0, meth=meth)
	foxOn_mean, foxOn_var = Sol_t[7], Sol_t[-1]
	# get the second moment from the variance
	# foxOn_s = foxOn_var + foxOn_mean**2
	return foxOn_mean, foxOn_var

def cell_pred_fn(k, times_arr, L, meth='BDF'):
	"""Function that gets the predictions for 6 concentrations each at 6 time points """
	nL_cons = len(L)
	nCons = int(nL_cons * len(times_arr))
	z0 = Get_Init_fn(k)
	# get the initial conditions here
	means_pred_arr = np.zeros(nCons)
	var_pred_arr = np.zeros(nCons)
	i = 0
	for igf in L:
		# ts = 0
		for t in times_arr:
			means_pred_arr[i], var_pred_arr[i] = FoxOn_preds_fn(k, igf, t, z0, meth=meth)
			i += 1
		# ts = t
	
	return means_pred_arr, var_pred_arr, z0[7]

# def get_percentiles
def get_prec_preds(Bounds_Perc_Dict, k, times_arr, L, meth='BDF'):
	"""" Gets the percentile prediction for one cell
	Inputs: Bounds_Perc_Dict: Dictionary of the bounds that produce specific percentiles
			k: parameter vector
			times_arr: times for solving the differential equation
			L: Ligand concentration for solving the differential equations

	"""
	
	means, var, init_point = cell_pred_fn(k, times_arr, L, meth=meth)
	alpha_arr = (means ** 2) / (var)
	scale_arr = var / means
	sumspp = np.zeros(len(alpha_arr))
	PredArr = np.array([])
	firstkey = list(Bounds_Perc_Dict.keys())[0]
	
	for key in Bounds_Perc_Dict:
		# print(key)
		bound = Bounds_Perc_Dict[key]
		if key == firstkey:
			
			spp = gamma.cdf(bound, a=alpha_arr, scale=scale_arr)  # specific percentile array
			sumspp += spp
		else:
			spp = gamma.cdf(bound, a=alpha_arr, scale=scale_arr) - sumspp
			sumspp += spp  # specific percentile array
		# print(spp)  #
		# the output means [percentile_L1_T1_bin1, percentile_L1_T2_bin1, percentile_L1_T3_bin1..... ...
		# .... percentile_L4_T6_bin9... percentile_L4_T7_bin9]
		#
		PredArr = np.concatenate((PredArr, spp), axis=0)
	return PredArr, means, var, init_point



