import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt

mdl_moments = np.loadtxt("/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Data/IGFR/Information_Values/Model_CeeMI_Step_Dose/mdl_moments_used.csv",delimiter=',')
expt_moments = np.loadtxt("/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Data/IGFR/Moments/experimental_moments_single_cell/multidose_moments.csv",delimiter=',')

mdl_mi = np.loadtxt("/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Data/IGFR/Information_Values/Model_CeeMI_Step_Dose/single_cell_mi_array.csv",delimiter=',')
expt_mi = np.loadtxt("/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Data/IGFR/Information_Values/Experimental_Step_Dose_MI_CSAR/Expt_CeeMI_Step_Dose.csv",delimiter=',')

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
fig,axs = plt.subplots(5,2,sharex = True)
colors = ["orange","green","blue","red"]
count = 0
for i in range(0,np.size(expt_mi),80):
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
for i in range(0,np.size(mdl_mi),10000):
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
fig.supxlabel("Response (Nuclear Foxo at Steady State)", x=0.5, y=0.05)
fig.supylabel("Probability")
plt.show()