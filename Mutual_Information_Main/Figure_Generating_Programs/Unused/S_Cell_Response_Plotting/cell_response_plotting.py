import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
import Mutual_Information_Main.functions.file_processing_functions.user_input_functions as uif

# Generates schematic showing cell state agnostic responses and single cell responses

plt.rcParams.update({'font.size': 15})

# _____ Load arguments BEGIN _____
# Specifies path to file which specifies user specified arguments, args file should be in same folder as this program
user_input_file_path = "args.txt"
user_arguments = uif.read_file_to_arg_dict(user_input_file_path)
# _____ Load arguments END _____

# _____ File path declarations BEGIN _____
# Specifies path to file which contains experimental moments
expt_moments_path = user_arguments["expt_moments_path"]

# Specifies path to file which contains moments from ode solution
mdl_moments_path = user_arguments["mdl_moments_path"]

# Specifies path to file which contains experimental MI values
expt_mi_path = user_arguments["expt_mi_path"]

# Specifies path to file which contains model MI values
mdl_mi_path = user_arguments["mdl_mi_path"]

# Specifies output path
output_path = user_arguments["output_path"]
# _____ File path declarations END _____

# _____ Loading files BEGIN _____
# Loads model moments
mdl_moments = np.loadtxt(mdl_moments_path,delimiter=',')
# Loads experimental moments
expt_moments = np.loadtxt(expt_moments_path,delimiter=',')
# Loads model mutual information values
mdl_mi = np.loadtxt(mdl_mi_path,delimiter=',')
# Loads experimental mutual information values
expt_mi = np.loadtxt(expt_mi_path,delimiter=',')
# _____ Loading files END _____



# _____ Figure generation BEGIN _____
info_sort_args = np.argsort(expt_mi)
expt_mi = expt_mi[info_sort_args]
expt_moments = expt_moments[info_sort_args,:]
expt_means = expt_moments[:, :4]
expt_vars = expt_moments[:, 4:] - expt_means ** 2
expt_scales = expt_vars / expt_means
expt_shapes = expt_means / expt_scales
min_resp = max(0, np.min(expt_means) - 2 * np.max(expt_vars ** .5))
max_resp = np.max(expt_means) + 2 * np.max(expt_vars ** .5)
x = np.linspace(min_resp, max_resp, 10000)
x = np.linspace(50, 700, 10000)
fig_dim = user_arguments["fig_dim"]
fig,axs = plt.subplots(20,2,sharex = True,figsize=(fig_dim[0], fig_dim[1]))
colors = user_arguments["colors"]
count = 0
for i in range(0,np.size(expt_mi),20):
	for j in range(4):
		x_sd = np.linspace(expt_means[i, j] - 3 * expt_vars[i, j] ** .5, expt_means[i, j] + 3 * expt_vars[i, j] ** .5, 1000)
		dist = st.gamma.pdf(x_sd, expt_shapes[i, j], scale=expt_scales[i, j])
		axs[count, 0].plot(x_sd, dist,color = colors[j])
	axs[0,0].set_xlim([50,700])
	axs[count, 0].spines["top"].set_visible(False)
	axs[count, 0].spines["right"].set_visible(False)
	axs[count, 0].spines["left"].set_visible(False)
	axs[count, 0].axes.yaxis.set_visible(False)
	#plt.title(f"I = {expt_mi[i]}\npercentile = {count/np.size(info_sort_args)*100}%")
	count += 1
axs[0, 0].set_title("Experiment")

info_sort_args = np.argsort(mdl_mi)
mdl_mi = mdl_mi[info_sort_args]
mdl_moments = mdl_moments[info_sort_args,:]
means = mdl_moments[:,:4]
vars = mdl_moments[:, 4:]-means**2
scales = vars/means
shapes = means/scales
min_resp = max(0,np.min(means) - 2*np.max(vars**.5))
max_resp = np.max(means) + 2*np.max(vars**.5)
#x = np.linspace(min_resp,max_resp,10000)
count = 0
for i in range(0,np.size(mdl_mi),2500):
	for j in range(4):
		x_sd = np.linspace(means[i, j]-3*vars[i,j]**.5,means[i, j]+3*vars[i,j]**.5,1000)
		dist = st.gamma.pdf(x_sd,shapes[i,j],scale = scales[i,j])
		axs[count, 1].plot(x_sd, dist,color = colors[j])
	#axs[count, 1].text(0.7,0.7,f"MI = {round(mdl_mi[i],3)}", horizontalalignment='center',verticalalignment='center', transform=axs[count, 1].transAxes)
	axs[count, 1].spines["top"].set_visible(False)
	axs[count, 1].spines["right"].set_visible(False)
	axs[count, 1].spines["left"].set_visible(False)
	axs[count, 1].axes.yaxis.set_visible(False)
	count += 1
axs[0, 1].set_title("Model")
fig.supxlabel("Response (Nuclear Foxo at Steady State)", x=0.5, y=0.02)
fig.supylabel("Probability")
plt.tight_layout()
plt.savefig(output_path)
plt.savefig(f"{output_path}.svg",format = 'svg')
plt.show()
# _____ Figure generation END _____
