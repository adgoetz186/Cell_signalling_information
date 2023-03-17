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
# Specifies path to file which specifies user specified arguments, args file should be in same folder as this program
user_input_file_path = "args.txt"
user_arguments = uif.read_file_to_arg_dict(user_input_file_path)
# _____ Load arguments END _____

# _____ File path declarations BEGIN _____
# Model multi-dose moments file path
mdl_md_moment_path = user_arguments["mdl_md_moment_path"]

# Model mutual information file path
mdl_mi_path = user_arguments["mdl_mi_path"]

# Path to single cell parameters
params = user_arguments["params"]
params = 10**np.loadtxt(params,delimiter=',')

# File path for figure image file
output_path = user_arguments["output_path"]
# _____ File path declarations END _____


# _____ Loading files BEGIN _____
# Loads Mutual information files, information taken from steady state responses
model_single_cell_mi = np.loadtxt(mdl_mi_path,delimiter=",")
model_nfoxo = np.loadtxt("/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Data/IGFR/pIGFR_pAkt_Ranges/nFoxo.csv",delimiter=",")
nfoxo_ranges = (model_nfoxo[:,4]+model_nfoxo[:,5])/2 - model_nfoxo[:,3]
print(nfoxo_ranges)
model_pakt = np.loadtxt("/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Data/IGFR/pIGFR_pAkt_Ranges/pAkt.csv",delimiter=",")
pakt_ranges = (model_pakt[:,4]+model_pakt[:,5])/2 - model_nfoxo[:,3]
print(pakt_ranges)
model_pigf = np.loadtxt("/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Data/IGFR/pIGFR_pAkt_Ranges/pIGF.csv",delimiter=",")
pigf_ranges = (model_pigf[:,4]+model_pigf[:,5])/2 - model_pigf[:,3]
print(pigf_ranges)
# Loads moment files

model_single_cell_moments = np.loadtxt(mdl_md_moment_path,delimiter=",")
print(st.pearsonr(nfoxo_ranges,model_single_cell_mi[:2000]))
print(st.pearsonr(pakt_ranges*params[:2000,8]/params[:2000,9],model_single_cell_mi[:2000]))
print(st.pearsonr(pigf_ranges*params[:2000,7]/params[:2000,6],model_single_cell_mi[:2000]))
input()
# _____ Loading files END _____


# _____ Data Processing BEGIN _____
# Generates response ranges by taking the difference in the mean responses
# from inputs of 0 pM to 125 pM at steady state (60 minutes after dose)

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

bin_count = user_arguments["bin_count"]
si_list = []
sensitivity_results = sa.bin_y_wrt_x(model_response_range, model_single_cell_mi,
                                     x_samp=response_range, y_samp=single_cell_mi, nbins=16,
                                     bin_center="center_of_bin_mass")
# sensitivity_results = sa.bin_y_wrt_x(np.arange(np.size(mdl_ifoxo))[0:400],model_single_cell_mi[0:400],x_samp = np.arange(np.size(single_cell_mean_initial_nuclear_foxo)),y_samp = single_cell_mi, nbins = 20, bin_center = "center_of_bin_mass")
# sensitivity_results = sa.bin_y_wrt_x(model_single_cell_mi[0:400],model_single_cell_mi[0:400],x_samp = single_cell_mi,y_samp = single_cell_mi, nbins = 20, bin_center = "center_of_bin_mass")
print(sensitivity_results)
var_mdl = sensitivity_results["total_pop_variance"]


var_exp_mdl = sensitivity_results["variance_of_expected_values_pop"]



exp_var_mdl = sensitivity_results["expected_value_of_variance_pop"]


# stores the expected value and locations of bins and the std of the samples within the bins
mdl_bin_center = sensitivity_results["mdl_bin_centers"]

mdl_information_bins = sensitivity_results['model_bin_heights']

mdl_information_bins_std_d = sensitivity_results['model_bin_std_d']

# _____ Data Processing END _____


# _____ Figure generation BEGIN _____
# sets figure colors
mdl_color = user_arguments["mdl_color"]

mdl_trend_color = user_arguments["mdl_trend_color"]


# initialize figure
fig, ax = plt.subplots(1)

# sets axis labels
ax.set_xlabel("Response Range")
ax.set_ylabel("Mutual Information (Bits)")

# allows kdeplot to have a standard legend label and visual
blue_patch = mpatches.Patch(color=mdl_color, label='Model',alpha = .6)

# Generates kernel density estimation plot for model data
sns.kdeplot(x=model_response_range, y=model_single_cell_mi,levels = [0.01,0.1,0.5,1],color = mdl_color,bw_adjust = 3,fill=True)

print(st.spearmanr(model_response_range,model_single_cell_mi))
# plots model trend line
lin_1 = plt.plot(mdl_bin_center,mdl_information_bins,label = "Model Trend", color = mdl_trend_color,linewidth = 3)


print(st.pearsonr(mdl_bin_center, mdl_information_bins))


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
plt.savefig(f"{output_path}.svg",format = 'svg')
plt.show()
# _____ Figure generation END _____
