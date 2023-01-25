import numpy as np

# Orders points based on x value, then performs binning of y values and generates important values for sensitivity
# analysis (primarily visual scatter-plot and analogs to Sobol indices).

def bin_y_wrt_x(x_pop,y_pop,x_samp = [],y_samp = [], nbins = 20, bin_center = "edge_center"):
	"""Orders points based on x value, then performs binning of y values and generates important values for sensitivity
	analysis
			Input:
				-x_pop (ndarray): array of points x values for the population values
				-y_pop (ndarray): array of points y values for the population values
				-x_samp (ndarray): array of points x values for sample values
				-y_pop (ndarray): array of points y values for sample values
				-y_pop (int): number of bins to use
				-bin_center (string): either 'edge_center' or 'center_of_bin_mass', specifies where to put the middle
					of the bin for visual plotting
			Output:
				dictionary of terms relevant for sensitivity analysis"""
	result_dict = {}
	
	# Finds bin edges while ensuring each bin gets a similar number of points
	percentile_edges = np.linspace(0, 100, nbins + 1)
	
	if x_samp != []:
		exp_bins = np.percentile(x_samp, percentile_edges)
	mdl_bins = np.percentile(x_pop, percentile_edges)
	
	
	# This adds a little to the upper bounds of the far right bin to avoid having a single point in the highest bin
	mdl_bins[-1] += 1
	if x_samp != []:
		exp_bins[-1] += 1
	
	# Finds which bin each point belongs to, binning is done based on response ranges
	mdl_rr_digi = np.digitize(x_pop, mdl_bins)
	if x_samp != []:
		exp_rr_digi = np.digitize(x_samp, exp_bins)
	
	if x_samp != []:
		# For the standard error bins the model using the same bins as used by the experimental data for comparison
		mdl_rr_exp_bin_digi = np.digitize(x_pop, exp_bins)
	
	# initializes numpy arrays to store bin data
	mdl_y_bins = np.zeros(nbins)
	mdl_y_bins_std_d = np.zeros(nbins)
	
	if x_samp != []:
		exp_y_bins = np.zeros(nbins)
		exp_y_bins_std_d = np.zeros(nbins)
		exp_y_bins_std_e = np.zeros(nbins)
		mdl_y_bins_std_e_ph = np.zeros(nbins)
		exp_mdl_y_bins_expt_def = np.zeros(nbins)
		total_mdl_data_expt_range = []
	
	for i in range(1, nbins + 1):
		# Calculates the binned average information values
		mdl_y_bins[i - 1] = np.mean(y_pop[np.argwhere(mdl_rr_digi == i)])
		print(np.size(y_pop[np.argwhere(mdl_rr_digi == i)]))
		if x_samp != []:
			exp_y_bins[i - 1] = np.mean(y_samp[np.argwhere(exp_rr_digi == i)])
		
		mdl_y_bins_std_d[i - 1] = np.std(y_pop[np.argwhere(mdl_rr_digi == i)], ddof=1)
		if x_samp != []:
			exp_y_bins_std_d[i - 1] = np.std(y_pop[np.argwhere(exp_rr_digi == i)], ddof=1)
			# Calculates the standard error for the sample data
			mdl_y_bins_std_e_ph[i - 1] = np.std(y_pop[np.argwhere(mdl_rr_exp_bin_digi == i)], ddof=1)
			print(exp_y_bins_std_e[i - 1])
			print(mdl_y_bins_std_e_ph[i - 1])
			print(np.size(y_samp[np.argwhere(exp_rr_digi == i)]))
			exp_y_bins_std_e[i - 1] = mdl_y_bins_std_e_ph[i - 1] / np.size(y_samp[np.argwhere(exp_rr_digi == i)]) ** .5
			exp_mdl_y_bins_expt_def[i-1] = np.mean(y_pop[np.argwhere(mdl_rr_exp_bin_digi == i)])
			total_mdl_data_expt_range.append(list(y_pop[np.argwhere(exp_rr_digi == i)]))
	if x_samp != []:
		result_dict["model_bin_std_d_expt_bins"] = mdl_y_bins_std_e_ph
		result_dict["model_bin_mean_expt_bins"] = exp_mdl_y_bins_expt_def
	result_dict["model_bin_heights"] = mdl_y_bins
	result_dict["model_bin_std_d"] = mdl_y_bins_std_d
	
	if x_samp != []:
		result_dict["expt_bin_heights"] = exp_y_bins
		result_dict["expt_bin_std_d"] = mdl_y_bins_std_d
		result_dict["expt_bin_std_e"] = exp_y_bins_std_e
		
	# Finds center of bins for plotting
	if bin_center == "edge_center":
		# This restores the true edge values
		mdl_bins[-1] -= 1
		if x_samp != []:
			exp_bins[-1] -= 1
		mdl_bin_center = (mdl_bins[1:] + mdl_bins[:-1]) / 2
		if x_samp != []:
			exp_bin_center = (exp_bins[1:] + exp_bins[:-1]) / 2
	
	if bin_center == "center_of_bin_mass":
		mdl_bin_center = np.zeros((nbins))
		if x_samp != []:
			exp_bin_center = np.zeros((nbins))
		for i in range(1, nbins + 1):
			# Calculates the binned average information values
			mdl_bin_center[i - 1] = np.mean(x_pop[np.argwhere(mdl_rr_digi == i)])
			if x_samp != []:
				exp_bin_center[i - 1] = np.mean(x_samp[np.argwhere(exp_rr_digi == i)])
	
	result_dict["mdl_bin_centers"] = mdl_bin_center
	if x_samp != []:
		result_dict["expt_bin_centers"] = exp_bin_center
	
	result_dict["total_pop_variance"] = np.var(y_pop,ddof=1)
	if x_samp != []:
		result_dict["total_samp_variance"] = np.var(y_samp, ddof=1)
		result_dict["total_pop_variance_sample_binned"] = np.var(np.array(total_mdl_data_expt_range), ddof=1)
		
	result_dict["variance_of_expected_values_pop"] = np.var(mdl_y_bins)
	if x_samp != []:
		result_dict["variance_of_expected_values_samp"] = np.var(exp_y_bins)
	
	if x_samp != []:
		result_dict["variance_of_expected_values_pop_expt_bin"] = np.var(exp_mdl_y_bins_expt_def)
	
	result_dict["expected_value_of_variance_pop"] = np.mean(mdl_y_bins_std_d**2)
	if x_samp != []:
		result_dict["expected_value_of_variance_samp"] = np.mean(exp_y_bins_std_d ** 2)
	
	return result_dict

			
	