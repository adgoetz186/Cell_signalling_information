import numpy as np
import json
import scipy.stats as st
import matplotlib.pyplot as plt

gillespie_results = np.loadtxt("/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Data/IGFR/Gillespie_Results/cell_0.csv",delimiter=',')
moment_results = np.loadtxt("/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Data/IGFR/Moments/Model_Moments/3_dose_response/moments_162_3_dose_full_time.csv",delimiter=',')
with open("/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Data/IGFR/Gillespie_Results/cell_header.txt",'r') as header:
	g_header = header.readlines()
g_header = [i.replace("'",'"') for i in g_header]
g_header = [json.loads(i) for i in g_header]
print(g_header)
print(gillespie_results)
print(g_header[:31])
dose_1 = gillespie_results[:,:31]
dose_2 = gillespie_results[:,31:62]
dose_3 = gillespie_results[:,62:]
time_list = list(set([i['time'] for i in g_header]))
time_list.sort()
print(time_list)


times = np.arange(0,90+3,3)*60
print(times)
print(np.size(times))
m_dose_1_mean = moment_results[:31]
m_dose_1_var = moment_results[93:124]-m_dose_1_mean**2
m_dose_1_scale = m_dose_1_var/m_dose_1_mean
m_dose_1_shape = m_dose_1_mean/m_dose_1_scale

m_dose_2_mean = moment_results[31:62]
m_dose_2_var = moment_results[124:155]-m_dose_2_mean**2
m_dose_2_scale = m_dose_2_var/m_dose_2_mean
m_dose_2_shape = m_dose_2_mean/m_dose_2_scale

m_dose_3_mean = moment_results[62:93]
m_dose_3_var = moment_results[155:186]-m_dose_3_mean**2
m_dose_3_scale = m_dose_3_var/m_dose_3_mean
m_dose_3_shape = m_dose_3_mean/m_dose_3_scale

fig, axs = plt.subplots(7,1,sharex=True)

time_indicies = [0,2,4,8,15,20,30]

for i in range(len(time_indicies)):
	
	x = np.linspace(np.min(dose_1[:,time_indicies[i]]),np.max(dose_1[:,time_indicies[i]]),10000)
	gamma_dist = st.gamma.pdf(x,m_dose_1_shape[time_indicies[i]],scale = m_dose_1_scale[time_indicies[i]])
	axs[i].plot(x,gamma_dist,color = (0,.7,1),zorder=10)
	
	x = np.linspace(np.min(dose_2[:, time_indicies[i]]), np.max(dose_2[:, time_indicies[i]]), 10000)
	gamma_dist = st.gamma.pdf(x, m_dose_2_shape[time_indicies[i]], scale=m_dose_2_scale[time_indicies[i]])
	axs[i].plot(x, gamma_dist,color = (1,.5,.5),zorder=11)
	
	
	x = np.linspace(np.min(dose_3[:, time_indicies[i]]), np.max(dose_3[:, time_indicies[i]]), 10000)
	gamma_dist = st.gamma.pdf(x, m_dose_3_shape[time_indicies[i]], scale=m_dose_3_scale[time_indicies[i]])
	axs[i].plot(x, gamma_dist,color = (.1,1,0),zorder=12)
	
	axs[i].hist(dose_1[:, time_indicies[i]], density='pdf',
	            bins=np.arange(np.min(dose_1[:, time_indicies[i]]), np.max(dose_1[:, time_indicies[i]]) + 1),
	            color="blue",zorder=10)
	axs[i].hist(dose_2[:, time_indicies[i]], density='pdf',
	            bins=np.arange(np.min(dose_2[:, time_indicies[i]]), np.max(dose_2[:, time_indicies[i]]) + 1),
	            color="orange",zorder=11)
	axs[i].hist(dose_3[:, time_indicies[i]], density='pdf',
	            bins=np.arange(np.min(dose_3[:, time_indicies[i]]), np.max(dose_3[:, time_indicies[i]]) + 1),
	            color="green",zorder=12)
	axs[i].spines["top"].set_visible(False)
	axs[i].spines["right"].set_visible(False)
	axs[i].spines["left"].set_visible(False)
	axs[i].axes.yaxis.set_visible(False)
	
plt.xlabel('Nuclear Foxo Levels')
plt.show()

#ax1.plot(time_list,np.mean(dose_1,axis=0))
#ax1.fill_between(time_list,np.mean(dose_1,axis=0)+np.std(dose_1,axis=0),np.mean(dose_1,axis=0)-np.std(dose_1,axis=0),alpha=0.4)
#ax1.plot(time_list,np.mean(dose_2,axis=0))
#ax1.fill_between(time_list,np.mean(dose_2,axis=0)+np.std(dose_2,axis=0),np.mean(dose_2,axis=0)-np.std(dose_2,axis=0),alpha=0.4)
#ax1.plot(time_list,np.mean(dose_3,axis=0))
#ax1.fill_between(time_list,np.mean(dose_3,axis=0)+np.std(dose_3,axis=0),np.mean(dose_3,axis=0)-np.std(dose_3,axis=0),alpha=0.4)
#ax1.spines["top"].set_visible(False)
#ax1.spines["right"].set_visible(False)


#ax2.plot(times,m_dose_1_mean)
#ax2.fill_between(times,m_dose_1_mean+m_dose_1_var**.5,m_dose_1_mean-m_dose_1_var**.5,alpha=0.4)
#ax2.plot(times,m_dose_2_mean)
#ax2.fill_between(times,m_dose_2_mean+m_dose_2_var**.5,m_dose_2_mean-m_dose_2_var**.5,alpha=0.4)
#ax2.plot(times,m_dose_3_mean)
#ax2.fill_between(times,m_dose_3_mean+m_dose_3_var**.5,m_dose_3_mean-m_dose_3_var**.5,alpha=0.4)
#ax2.spines["top"].set_visible(False)
#ax2.spines["right"].set_visible(False)

#plt.show()