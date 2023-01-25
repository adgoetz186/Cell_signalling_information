import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
gillespie_results = np.loadtxt("/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Data/EGFR/Gillespie_Results/cell_0.csv",delimiter=',')
mdl_means = np.loadtxt("/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Data/EGFR/means_and_variances/noise_factor_1/mu_226.csv",delimiter=',')[0,:]
mdl_vars = np.loadtxt("/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Data/EGFR/means_and_variances/noise_factor_1/var_226.csv",delimiter=',')[0,:]
with open("/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Data/EGFR/Gillespie_Results/cell_header.txt",'r') as header:
	g_header = header.readlines()

fig,axs = plt.subplots(10,1,sharex=True)
for i in range(np.size(mdl_means)):
	axs[i].hist(gillespie_results[:,i]-np.mean(gillespie_results[:,i]),density='pdf',bins=20)
	axs[i].spines["top"].set_visible(False)
	axs[i].spines["right"].set_visible(False)
	axs[i].spines["left"].set_visible(False)
	axs[i].axes.yaxis.set_visible(False)

mdl_scale = mdl_vars/mdl_means
mdl_shape = mdl_means/mdl_scale
mins = np.min(gillespie_results,axis=0)
maxs = np.max(gillespie_results,axis=0)


for i in range(np.size(mdl_means)):
	shift = np.mean(gillespie_results[:,i])
	x = np.linspace(mins[i]-shift, maxs[i]-shift, 10000)
	pdf = st.gamma.pdf(x,mdl_shape[i],scale = mdl_scale[i],loc = -shift)
	axs[i].plot(x,pdf)
plt.xlabel("sEGFR at Steady State")
plt.xticks(ticks = [-1000,0,1000],labels = ["$\mu$-1000","$\mu$","$\mu$+1000"])
plt.show()