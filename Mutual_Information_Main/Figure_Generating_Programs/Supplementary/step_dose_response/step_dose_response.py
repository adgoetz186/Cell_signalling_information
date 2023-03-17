import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import os

# sets font for figure text
plt.rcParams.update({'font.size': 15})

# _____ Setting the CWD to be Mutual_Information_Main BEGIN _____
# Cell_signaling_information path here
path_to_CSI = ""
if path_to_CSI == "":
	try:
		# Obtains the location of the Cell_signaling_information folder if it is in the cwd parents
		path_to_CSI = Path.cwd().parents[[Path.cwd().parents[i].parts[-1] for i in range(len(Path.cwd().parts) - 1)].index(
			"Cell_signalling_information")]
	except ValueError:
		print("Cell_signalling_information not found in cwd parents, trying sys.path")
		try:
			# Obtains the location of the Cell_signaling_information folder if it is in sys.path
			path_to_CSI = Path(sys.path[[Path(i).parts[-1] for i in sys.path].index("Cell_signalling_information")])
		except ValueError:
			print("Cell_signalling_information not found in sys.path "
			      "consult 'Errors with setting working directory' in README")
else:
	path_to_CSI = Path(path_to_CSI)
path_to_MIM = path_to_CSI/"Mutual_Information_Main"
os.chdir(path_to_MIM)
# _____ Setting the CWD to be Mutual_Information_Main END _____

multidose_location = Path("Data/IGFR/experimental_data/adjusted_multidose/multi_dose.npy")
output_location = Path("Figures/Supplementary/step_dose_response/")
output_name = "step_dose"
multidose = np.load(multidose_location)
times = multidose[0,:]
doses = multidose[1:,:]
average_dose = np.average(doses,axis = 0)
fig,(ax1,ax2) = plt.subplots(1,2,figsize=(12, 4))

ax1.plot(times,np.transpose(doses[1:,:]),color = 'grey',alpha = 0.3)
ax1.plot(times,doses[0],color = 'grey',alpha = 0.3,label = "Single Cell")
#ax.plot(times,np.transpose(doses),color = 'orange',alpha = 0.3)
ax1.plot(times,average_dose,color = 'k', label = 'Average')

ymin = 50
ymax = 850
fct = times[0:11]
sct = times[30:41]
tct = times[60:71]
foct = times[90:101]
ax1.fill_between(fct, np.ones_like(fct)*(ymin-10), np.ones_like(fct)*(ymax+10), color='#FAA0A0', alpha=.7,zorder=3)
ax1.fill_between(sct, np.ones_like(sct)*(ymin-10), np.ones_like(sct)*(ymax+10), color='#C1E1C1', alpha=.7,zorder=3)
ax1.fill_between(tct, np.ones_like(tct)*(ymin-10), np.ones_like(tct)*(ymax+10), color='#AEC6CF', alpha=.7,zorder=3)
ax1.fill_between(foct, np.ones_like(foct)*(ymin-10), np.ones_like(foct)*(ymax+10), color='#C3B1E1', alpha=.7,zorder=3)
#ax.fill_between(times,np.percentile(doses,25,axis = 0),np.percentile(doses,75,axis = 0),alpha = 0.5,color = 'blue', label = 'spread')
ax1.set_ylim([ymin,ymax])
ax1.set_ylabel("Nuclear FoxO")
ax1.set_xlabel("Time (min)")
ax1.legend(frameon=False,ncol=2,loc="upper center")
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

ax2.plot(times,doses[0],color = 'black',label = "Single Cell")

ymin = 175
ymax = 500
fct = times[0:11]
sct = times[30:41]
tct = times[60:71]
foct = times[90:101]
ax2.fill_between(fct, np.ones_like(fct)*(ymin-10), np.ones_like(fct)*(ymax+10), color='#FAA0A0', alpha=.7,zorder=3)
ax2.fill_between(sct, np.ones_like(sct)*(ymin-10), np.ones_like(sct)*(ymax+10), color='#C1E1C1', alpha=.7,zorder=3)
ax2.fill_between(tct, np.ones_like(tct)*(ymin-10), np.ones_like(tct)*(ymax+10), color='#AEC6CF', alpha=.7,zorder=3)
ax2.fill_between(foct, np.ones_like(foct)*(ymin-10), np.ones_like(foct)*(ymax+10), color='#C3B1E1', alpha=.7,zorder=3)
#ax.fill_between(times,np.percentile(doses,25,axis = 0),np.percentile(doses,75,axis = 0),alpha = 0.5,color = 'blue', label = 'spread')
ax2.set_ylim([ymin,ymax])
ax2.set_ylabel("Nuclear FoxO")
ax2.set_xlabel("Time (min)")
ax2.legend(frameon=False,ncol=2,loc="upper center")
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig(output_location/output_name)
plt.savefig(output_location/(output_name+".svg"),format = 'svg')
plt.show()