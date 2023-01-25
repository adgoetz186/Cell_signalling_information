import numpy as np
import seaborn as sns
import random
import scipy.stats as st
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import Mutual_Information_Final_Version.functions.data_processing_functions.sensitivity_analysis as sa
import Mutual_Information_Final_Version.functions.file_processing_functions.user_input_functions as uif
plt.rcParams.update({'font.size': 15})

# Generates scatter plot for MI vs initial nuclear foxo

# _____ Load arguments BEGIN _____
# Specifies path to file which specifies user specified arguments
user_input_file_path = "/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Figure_Generating_Programs/Figure_3/E/args.txt"
user_arguments = uif.read_file_to_arg_dict(user_input_file_path)
# _____ Load arguments END _____

# _____ File path declarations BEGIN _____
# Experimental multi-dose trajectories file path
expt_md_moment_path = user_arguments["expt_md_moment_path"]

# Experimental multi-dose trajectories file path
params = "/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Data/IGFR/Information_Values/Model_CeeMI_Step_Dose/mdl_params_used.csv"
params = 10**np.loadtxt(params,delimiter=',')
# Model multi-dose moments file path
mdl_md_moment_path = user_arguments["mdl_md_moment_path"]

# Experimental mutual information file path
expt_mi_path = user_arguments["expt_mi_path"]

# Model mutual information file path
mdl_mi_path = user_arguments["mdl_mi_path"]

# File path for figure image file
output_path = user_arguments["output_path"]
# _____ File path declarations END _____


# _____ Loading files BEGIN _____
# Loads Mutual information files, information taken from steady state responses
single_cell_mi = np.loadtxt(expt_mi_path, delimiter=",")
model_single_cell_mi = np.loadtxt(mdl_mi_path,delimiter=",")

# Loads experimental trajectory file
single_cell_nuclear_foxo = np.loadtxt(expt_md_moment_path,delimiter=",")

# Loads model moment file
model_single_cell_moments = np.loadtxt(mdl_md_moment_path,delimiter=",")
# _____ Loading files END _____


# _____ Data Processing BEGIN _____
# Finds each cells initial nuclear foxo averaged over all recorded pre dose measurements for each cell
single_cell_mean_initial_nuclear_foxo = np.average(single_cell_nuclear_foxo[1:,:11],axis = 1)

# Model initial nuclear foxo
mdl_ifoxo = model_single_cell_moments[:,0]
nfoxo = np.loadtxt("/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Data/IGFR/Gillespie_Results/initial_foxo.csv")
print(nfoxo)

for i in range(400):
	print(i)
	print(mdl_ifoxo[i])
print(st.spearmanr(nfoxo,model_single_cell_mi[0:400]))
print(st.spearmanr(mdl_ifoxo[0:400],model_single_cell_mi[0:400]))

var_val_ar = 10 ** np.linspace(-3, 3, 20)
pc_list = []
tot_val = 0
div = 0
#for var_val in var_val_ar:
#	val = 0
#	for i in range(1000):
#		print(var_val, i)
#		index_model_cells = np.arange(np.shape(model_single_cell_moments)[0])
#		np.random.shuffle(index_model_cells)
#		index_model_cells = index_model_cells[0:400]
#		index_subsamples_range = []
#		for imc in range(len(index_model_cells)):
#			low_var = model_single_cell_moments[index_model_cells[imc], 4] - model_single_cell_moments[
#				index_model_cells[imc], 0] ** 2
#			print(low_var,model_single_cell_moments[index_model_cells[imc], 0])
#			low_var /= var_val
#
#			low_scale = low_var / model_single_cell_moments[index_model_cells[imc], 0]
#			low_shape = model_single_cell_moments[index_model_cells[imc], 0] / low_scale
#			low_response = st.gamma.rvs(low_shape, scale=low_scale)
#
#			index_subsamples_range.append(low_response)
#		tot_val+=st.spearmanr(model_single_cell_moments[index_model_cells, 0], model_single_cell_mi[index_model_cells])[0]
#		div +=1
#		model_single_cell_mi_samp = model_single_cell_mi[index_model_cells]
#		val += st.spearmanr(index_subsamples_range, model_single_cell_mi_samp)[0]
#	pc_list.append(val / 1000)
#plt.plot(var_val_ar, pc_list, label="model sample with sampled response", color='green')
print(st.pearsonr(single_cell_mean_initial_nuclear_foxo, single_cell_mi))
print(st.pearsonr(mdl_ifoxo[0:400], model_single_cell_mi[0:400]))
#plt.hlines(st.spearmanr(single_cell_mean_initial_nuclear_foxo, single_cell_mi)[0], np.min(var_val_ar),
#           np.max(var_val_ar),
#           label="experimental", color='red')
#plt.hlines(tot_val/div, np.min(var_val_ar),
#           np.max(var_val_ar), label="model sample with mean response", color='blue')
#plt.legend()
#plt.xlabel("sampled variance = stochastic variance / n")
#plt.ylabel("Spearman cor")
#plt.xscale("log")
#plt.show()


bin_count = user_arguments["bin_count"]
K_dis = params[:,3]/params[:,2]
K_ip = params[:,4]/params[:,5]
K_ap = params[:,7]/params[:,6]
K_fp = params[:,8]/params[:,9]
K_fio = params[:,11]/params[:,10]
sensitivity_results_dis = sa.bin_y_wrt_x(np.log10(K_dis),model_single_cell_mi,nbins = 16, bin_center = "center_of_bin_mass")
sensitivity_results_ip = sa.bin_y_wrt_x(np.log10(K_ip),model_single_cell_mi,nbins = 16, bin_center = "center_of_bin_mass")
sensitivity_results_ap = sa.bin_y_wrt_x(np.log10(K_ap),model_single_cell_mi,nbins = 16, bin_center = "center_of_bin_mass")
sensitivity_results_fp = sa.bin_y_wrt_x(np.log10(K_fp),model_single_cell_mi,nbins = 16, bin_center = "center_of_bin_mass")
sensitivity_results_fio = sa.bin_y_wrt_x(np.log10(K_fio),model_single_cell_mi,nbins = 16, bin_center = "center_of_bin_mass")

print(bin_count)
# Finds bin edges while ensuring each bin gets a similar number of points
percentile_edges = np.linspace(0, 100, bin_count + 1)
mdl_bins = np.percentile(mdl_ifoxo, percentile_edges)
exp_bins = np.percentile(single_cell_mean_initial_nuclear_foxo, percentile_edges)

# Finds center of bins for plotting
mdl_bin_center = (mdl_bins[1:] + mdl_bins[:-1]) / 2
exp_bin_center = (exp_bins[1:] + exp_bins[:-1]) / 2

# This adds a little to the upper bounds of the far right bin to avoid having a single point in the highest bin
mdl_bins[-1] += 1
exp_bins[-1] += 1

# Finds which bin each point belongs to, binning is done based on initial nuclear foxo
mdl_ifoxo_digi = np.digitize(mdl_ifoxo, mdl_bins)
exp_ifoxo_digi = np.digitize(single_cell_mean_initial_nuclear_foxo, exp_bins)

# For the standard error bins the model using the same bins as used by the experimental data for comparison
mdl_ifoxo_exp_bin_digi = np.digitize(mdl_ifoxo, exp_bins)

# initializes numpy arrays to store bin data
mdl_information_bins = np.zeros(bin_count)
exp_information_bins = np.zeros(bin_count)
mdl_information_bins_std = np.zeros(bin_count)
exp_information_bins_std = np.zeros(bin_count)

for i in range(1, bin_count + 1):
	# Calculates the binned average information values
	mdl_information_bins[i - 1] = np.average(model_single_cell_mi[np.argwhere(mdl_ifoxo_digi == i)])
	exp_information_bins[i - 1] = np.average(single_cell_mi[np.argwhere(exp_ifoxo_digi == i)])
	
	# Calculates the standard error for the sample data
	mdl_information_bins_std[i - 1] = np.std(model_single_cell_mi[np.argwhere(mdl_ifoxo_exp_bin_digi == i)])
	exp_information_bins_std[i - 1] = mdl_information_bins_std[i - 1] / np.size(
		single_cell_mi[np.argwhere(exp_ifoxo_digi == i)]) ** .5
# _____ Data Processing END _____


# _____ Figure generation BEGIN _____
# sets the figure colors
mdl_color = user_arguments["mdl_color"]
expt_color = user_arguments["expt_color"]
mdl_trend_color = user_arguments["mdl_trend_color"]
expt_trend_color = user_arguments["expt_trend_color"]

fig,ax = plt.subplots(1)
sns.kdeplot(x=K_fio, y=model_single_cell_mi, levels=[0.01, 0.1, 0.5, 1], color=mdl_color,bw_adjust=3, fill=True, alpha=1,log_scale=(10,False))
#plt.errorbar(10**sensitivity_results_fio["mdl_bin_centers"],sensitivity_results_fio["model_bin_heights"],yerr=sensitivity_results_fio["model_bin_std_d"],color = mdl_trend_color,linewidth = 3)
ax.plot(10**sensitivity_results_fio["mdl_bin_centers"],sensitivity_results_fio["model_bin_heights"],color = mdl_trend_color,linewidth = 3)
ax.set_ylim([1,2.01])
ax.set_xlim([.014,3.1])
ax.set_xscale("log")
ax.set_ylabel("Mutual Information (bits)")
ax.set_xlabel("Initial Cytoplasm/Nuclear Foxo Ratio")
pearsonr = st.pearsonr(np.log10(K_fio),model_single_cell_mi)
#ax.set_title(f"$\\rho$ = {np.round(pearsonr[0],4)}, log(pval) = {np.round(np.log10(pearsonr[1]),2)}")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.show()


fig,ax = plt.subplots(1)
sns.kdeplot(x=K_dis, y=model_single_cell_mi, levels=[0.01, 0.1, 0.5, 1], color=mdl_color,bw_adjust=3, fill=True, alpha=1,log_scale=(10,False))
#plt.scatter(np.log10(K_dis),model_single_cell_mi,s = 1,color = 'red')
#plt.errorbar(10**sensitivity_results_dis["mdl_bin_centers"],sensitivity_results_dis["model_bin_heights"],yerr=sensitivity_results_dis["model_bin_std_d"],color = mdl_trend_color,linewidth = 3)
ax.plot(10**sensitivity_results_dis["mdl_bin_centers"],sensitivity_results_dis["model_bin_heights"],color = mdl_trend_color,linewidth = 3)
ax.set_ylim([1,2.01])
ax.set_xlim([3,750])
#plt.xscale("log")
ax.set_ylabel("Mutual Information (bits)")
ax.set_xlabel("IGFR Dissociation Constant")
pearsonr = st.pearsonr(np.log10(K_dis),model_single_cell_mi)
#plt.title(f"$\\rho$ = {np.round(pearsonr[0],4)}, log(pval) = {np.round(np.log10(pearsonr[1]),2)}")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.show()


fig,ax = plt.subplots(1)
sns.kdeplot(x=K_ip, y=model_single_cell_mi, levels=[0.01, 0.1, 0.5, 1], color=mdl_color,bw_adjust=3, fill=True, alpha=1,log_scale=(10,False))
#plt.errorbar(10**sensitivity_results_ip["mdl_bin_centers"],sensitivity_results_ip["model_bin_heights"],yerr=sensitivity_results_ip["model_bin_std_d"],color = mdl_trend_color,linewidth = 3)
ax.plot(10**sensitivity_results_ip["mdl_bin_centers"],sensitivity_results_ip["model_bin_heights"],color = mdl_trend_color,linewidth = 3)
ax.set_ylim([1,2.01])
ax.set_xlim([.08,16])
ax.set_xscale("log")
ax.set_ylabel("Mutual Information (bits)")
ax.set_xlabel("$K_{pIGF}$")
pearsonr = st.pearsonr(np.log10(K_ip),model_single_cell_mi)
#plt.title(f"$\\rho$ = {np.round(pearsonr[0],4)}, log(pval) = {np.round(np.log10(pearsonr[1]),2)}")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.show()


fig,ax = plt.subplots(1)
sns.kdeplot(x=K_ap, y=model_single_cell_mi, levels=[0.01, 0.1, 0.5, 1], color=mdl_color,bw_adjust=3, fill=True, alpha=1,log_scale=(10,False))
#plt.errorbar(10**sensitivity_results_ap["mdl_bin_centers"],sensitivity_results_ap["model_bin_heights"],yerr=sensitivity_results_ap["model_bin_std_d"],color = mdl_trend_color,linewidth = 3)
ax.plot(10**sensitivity_results_ap["mdl_bin_centers"],sensitivity_results_ap["model_bin_heights"],color = mdl_trend_color,linewidth = 3)
ax.set_ylim([1,2.01])
ax.set_xlim([5e-5,.01])
ax.set_xscale("log")
ax.set_ylabel("Mutual Information (bits)")
ax.set_xlabel("$K_{pAkt}$")
pearsonr = st.pearsonr(np.log10(K_ap),model_single_cell_mi)
#plt.title(f"$\\rho$ = {np.round(pearsonr[0],4)}, log(pval) = {np.round(np.log10(pearsonr[1]),2)}")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.show()


fig,ax = plt.subplots(1)
sns.kdeplot(x=K_fp, y=model_single_cell_mi, levels=[0.01, 0.1, 0.5, 1], color=mdl_color,bw_adjust=3, fill=True, alpha=1,log_scale=(10,False))
#plt.errorbar(10**sensitivity_results_fp["mdl_bin_centers"],sensitivity_results_fp["model_bin_heights"],yerr=sensitivity_results_fp["model_bin_std_d"],color = mdl_trend_color,linewidth = 3)
ax.plot(10**sensitivity_results_fp["mdl_bin_centers"],sensitivity_results_fp["model_bin_heights"],color = mdl_trend_color,linewidth = 3)
ax.set_ylim([1,2.01])
ax.set_xlim([.0004,.08])
ax.set_xscale("log")
ax.set_ylabel("Mutual Information (bits)")
ax.set_xlabel("$k_{pFoxO}$")
pearsonr = st.pearsonr(np.log10(K_fp),model_single_cell_mi)
#plt.title(f"$\\rho$ = {np.round(pearsonr[0],4)}, log(pval) = {np.round(np.log10(pearsonr[1]),2)}")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.show()




# initialize figure
fig, ax = plt.subplots(1)

# sets axis labels
ax.set_xlabel("Initial nuclear Foxo")
ax.set_ylabel("Mutual Information (Bits)")

# allows kdeplot to have a standard legend label and visual
blue_patch = mpatches.Patch(color=mdl_color, label='Model', alpha=.6)

# Generates kernel density estimation plot for model data
sns.kdeplot(x=mdl_ifoxo, y=model_single_cell_mi, levels=[0.01, 0.1, 0.5, 1], color=mdl_color,bw_adjust=3, fill=True, alpha=1)

# Generates scatter plot for experimental data
val_2 = plt.scatter(single_cell_mean_initial_nuclear_foxo,single_cell_mi,label = "Experiment",color = expt_color)

print(st.spearmanr(single_cell_mean_initial_nuclear_foxo,single_cell_mi))
print(st.spearmanr(mdl_ifoxo,model_single_cell_mi))
# removes plot spines
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# plots model trend line
lin_1 = plt.plot(mdl_bin_center,mdl_information_bins,label = "Model Trend", color = mdl_trend_color,linewidth = 3)

# plots experimental trend line
lin_2 = plt.errorbar(exp_bin_center,exp_information_bins,yerr=exp_information_bins_std,label = "Expt Trend", color = expt_trend_color,linewidth = 3)

print(st.pearsonr(mdl_bin_center, mdl_information_bins))
print(st.pearsonr(exp_bin_center, exp_information_bins))

# generates legend
ax.legend(handles=[blue_patch,val_2,lin_1[0],lin_2], loc="lower right",frameon=False)

# tweaks axis properties for visual effect
ax.set_ylim([1,2.02])
ax.set_xlim([290,780])

plt.tight_layout()
plt.savefig(output_path)
plt.show()
# _____ Figure generation END _____
