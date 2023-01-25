import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import Mutual_Information_Final_Version.functions.data_processing_functions.sensitivity_analysis as sa
import random
import scipy.stats as st
import Mutual_Information_Final_Version.functions.file_processing_functions.user_input_functions as uif
# Generates scatter plot for MI vs experimental response range

# sets font for figure text
plt.rcParams.update({'font.size': 15})

# _____ Load arguments BEGIN _____
# Specifies path to file which specifies user specified arguments
user_input_file_path = "/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Figure_Generating_Programs/Figure_3/D/args.txt"
user_arguments = uif.read_file_to_arg_dict(user_input_file_path)
# _____ Load arguments END _____

# _____ File path declarations BEGIN _____
# Experimental multi-dose moments file path
expt_md_moment_path = user_arguments["expt_md_moment_path"]

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

# Loads moment files
single_cell_moments = np.loadtxt(expt_md_moment_path,delimiter=",")
model_single_cell_moments = np.loadtxt(mdl_md_moment_path,delimiter=",")
# _____ Loading files END _____


# _____ Data Processing BEGIN _____
# Generates response ranges by taking the difference in the mean responses
# from inputs of 0 pM to 125 pM at steady state (60 minutes after dose)
response_range = single_cell_moments[:,0]-single_cell_moments[:,3]
#response_range = np.max(single_cell_moments[:,:4],axis = 1)-np.min(single_cell_moments[:,:4],axis = 1)
model_response_range = model_single_cell_moments[:,0]-model_single_cell_moments[:,3]

var_val_ar = 10**np.linspace(-3,3,20)
pc_list = []
#for var_val in var_val_ar:
#	val = 0
#	for i in range(1000):
#		print(var_val,i)
#		index_model_cells = np.arange(np.shape(model_single_cell_moments)[0])
#		np.random.shuffle(index_model_cells)
#		index_model_cells = index_model_cells[0:400]
#		index_subsamples_range = []
#		for imc in range(len(index_model_cells)):
#			low_var = model_single_cell_moments[index_model_cells[imc], 4] - model_single_cell_moments[index_model_cells[imc], 0]**2
#			low_var/=var_val
#			low_scale = low_var/model_single_cell_moments[index_model_cells[imc], 0]
#			low_shape = model_single_cell_moments[index_model_cells[imc], 0]/low_scale
#			low_response = st.gamma.rvs(low_shape,scale = low_scale)
#
#			high_var = model_single_cell_moments[index_model_cells[imc], 7] - model_single_cell_moments[index_model_cells[imc], 3] ** 2
#			high_var /= var_val
#			high_scale = high_var / model_single_cell_moments[index_model_cells[imc], 3]
#			high_shape = model_single_cell_moments[index_model_cells[imc], 3] / high_scale
#			high_response = st.gamma.rvs(high_shape, scale=high_scale)
#
#			index_subsamples_range.append(low_response-high_response)
#
#		model_single_cell_mi_samp = model_single_cell_mi[index_model_cells]
#		val += st.spearmanr(index_subsamples_range,model_single_cell_mi_samp)[0]
#	pc_list.append(val/1000)
#plt.plot(var_val_ar,pc_list,label ="model sample with sampled response",color = 'green')
print(st.pearsonr(response_range,single_cell_mi))
print(st.pearsonr(model_response_range[0:400],model_single_cell_mi[0:400]))
#plt.hlines(st.spearmanr(response_range,single_cell_mi)[0],np.min(var_val_ar),np.max(var_val_ar),label = "experimental",color = 'red')
#plt.hlines(st.spearmanr(model_response_range[0:400],model_single_cell_mi[0:400])[0],np.min(var_val_ar),np.max(var_val_ar),label = "model sample with mean response",color = 'blue')
#plt.legend()
#plt.xlabel("sampled variance = stochastic variance / n")
#plt.ylabel("Spearman cor")
#plt.xscale("log")
#plt.show()
#input()
plt.scatter(response_range, single_cell_mi)
count = 0
for i in range(np.shape(single_cell_moments)[0]):
    
    print(i)
    print(single_cell_moments[i])
    diff = single_cell_moments[i,:3] - single_cell_moments[i,1:4]
    print(diff)
    if np.sum(diff) != np.sum(np.abs(diff)):
        count +=1
        plt.scatter(response_range[i],single_cell_mi[i],color = "red",s = 50)
        print(single_cell_moments[i,0],single_cell_moments[i,1],single_cell_moments[i,2],single_cell_moments[i,3])
print(count)
plt.ylabel("Mutual Information")
plt.xlabel("Experimental response range")
plt.show()
bin_count = user_arguments["bin_count"]
si_list = []
sensitivity_results = sa.bin_y_wrt_x(model_response_range, model_single_cell_mi,
                                     x_samp=response_range, y_samp=single_cell_mi, nbins=16,
                                     bin_center="center_of_bin_mass")
# sensitivity_results = sa.bin_y_wrt_x(np.arange(np.size(mdl_ifoxo))[0:400],model_single_cell_mi[0:400],x_samp = np.arange(np.size(single_cell_mean_initial_nuclear_foxo)),y_samp = single_cell_mi, nbins = 20, bin_center = "center_of_bin_mass")
# sensitivity_results = sa.bin_y_wrt_x(model_single_cell_mi[0:400],model_single_cell_mi[0:400],x_samp = single_cell_mi,y_samp = single_cell_mi, nbins = 20, bin_center = "center_of_bin_mass")
print(sensitivity_results)
var_mdl = sensitivity_results["total_pop_variance"]
var_expt = sensitivity_results["total_samp_variance"]

var_exp_mdl = sensitivity_results["variance_of_expected_values_pop"]
var_exp_expt = sensitivity_results["variance_of_expected_values_samp"]
var_exp_mdl_expt_bin = sensitivity_results["variance_of_expected_values_pop_expt_bin"]

exp_var_mdl = sensitivity_results["expected_value_of_variance_pop"]
exp_var_expt = sensitivity_results["expected_value_of_variance_samp"]
print(var_exp_mdl_expt_bin/sensitivity_results["total_pop_variance_sample_binned"])
print(var_exp_mdl / var_mdl, var_exp_expt / var_expt)
print(exp_var_mdl / var_mdl, exp_var_expt / var_expt)
print(var_exp_mdl / exp_var_mdl, var_exp_expt / exp_var_expt)
input()
# stores the expected value and locations of bins and the std of the samples within the bins
mdl_bin_center = sensitivity_results["mdl_bin_centers"]
exp_bin_center = sensitivity_results["expt_bin_centers"]
mdl_information_bins = sensitivity_results['model_bin_heights']
exp_information_bins = sensitivity_results['expt_bin_heights']
mdl_information_bins_std_d = sensitivity_results['model_bin_std_d']
exp_information_bins_std_e = sensitivity_results['expt_bin_std_e']
# _____ Data Processing END _____


# _____ Figure generation BEGIN _____
# sets figure colors
mdl_color = user_arguments["mdl_color"]
expt_color = user_arguments["expt_color"]
mdl_trend_color = user_arguments["mdl_trend_color"]
expt_trend_color = user_arguments["expt_trend_color"]

# initialize figure
fig, ax = plt.subplots(1)

# sets axis labels
ax.set_xlabel("Response Range")
ax.set_ylabel("Mutual Information (Bits)")

# allows kdeplot to have a standard legend label and visual
blue_patch = mpatches.Patch(color=mdl_color, label='Model',alpha = .6)

# Generates kernel density estimation plot for model data
sns.kdeplot(x=model_response_range, y=model_single_cell_mi,levels = [0.01,0.1,0.5,1],color = mdl_color,bw_adjust = 3,fill=True)


# Generates scatter plot for experimental data
val = plt.scatter(response_range,single_cell_mi,label = "Experiment",color = expt_color)
print(st.spearmanr(response_range,single_cell_mi))
print(st.spearmanr(model_response_range,model_single_cell_mi))
# plots model trend line
lin_1 = plt.plot(mdl_bin_center,mdl_information_bins,label = "Model Trend", color = mdl_trend_color,linewidth = 3)

# plots experimental trend line
lin_2 = plt.errorbar(exp_bin_center,exp_information_bins,yerr=exp_information_bins_std_e,label = "Expt Trend", color = expt_trend_color,linewidth = 3)

print(st.pearsonr(mdl_bin_center, mdl_information_bins))
print(st.pearsonr(exp_bin_center, exp_information_bins))

# generates legend
#ax.legend(handles=[blue_patch,lin_1[0],val,lin_2],loc="upper center",frameon=False,ncol = 2, bbox_to_anchor=(.5,1.25))

# removes plot spines
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# tweaks axis properties for visual effect
ax.set_ylim([1,2.02])
ax.set_yticks([1,1.2,1.4,1.6,1.8,2])
ax.set_xticks([100,200,300,400,500])
ax.set_xlim([90,480])

img = plt.imread("/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Figures/Figure_3/D/range.png")
# above legend
#newax = fig.add_axes([0.65,0.45,0.3,0.3], anchor='E', zorder=1)
# left of legend
#newax = fig.add_axes([0.3,0.15,0.3,0.3], anchor='E', zorder=1)
# at legend
newax = fig.add_axes([0.55,0.18,0.35,0.35], anchor='E', zorder=1)
newax.imshow(img)
newax.axis('off')

# generates and saves figure in the tight layout
plt.tight_layout()
plt.savefig(output_path)
plt.show()
# _____ Figure generation END _____
