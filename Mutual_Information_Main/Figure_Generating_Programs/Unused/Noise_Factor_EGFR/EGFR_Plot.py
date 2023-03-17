import numpy as np
import matplotlib.pyplot as plt
import Mutual_Information_Final_Version.functions.file_processing_functions.user_input_functions as uif

# Generates plot of EGFR signaling performance with difference levels of noise and MI of cell state agnostic response

# sets font for figure text
plt.rcParams.update({'font.size': 15})

# _____ Load arguments BEGIN _____
# Specifies path to file which specifies user specified arguments, args file should be in same folder as this program
user_input_file_path = "args.txt"
user_arguments = uif.read_file_to_arg_dict(user_input_file_path)
# _____ Load arguments END _____

# _____ File path declarations BEGIN _____
# location of experimental CSAR MI
CSAR_mi = user_arguments["CSAR_mi"]
# location of nf 1 single cell MI values
nf1_sc_mi = user_arguments["nf1_sc_mi"]
# location of nf 2 single cell MI values
nf2_sc_mi = user_arguments["nf2_sc_mi"]
# location of nf 4 single cell MI values
nf4_sc_mi = user_arguments["nf4_sc_mi"]
# location of nf 8 single cell MI values
nf8_sc_mi = user_arguments["nf8_sc_mi"]
# location of nf 16 single cell MI values
nf16_sc_mi = user_arguments["nf16_sc_mi"]
# location of nf 32 single cell MI values
nf32_sc_mi = user_arguments["nf32_sc_mi"]
# location of output image
output_dir = user_arguments["output_dir"]
# _____ File path declarations END _____

# _____ Loading files BEGIN _____
# loads the CSAR mutual information performance value
CSAR_mi_array = np.loadtxt(CSAR_mi)

# loads single cell mutual information performance values for various noise factors
nf1_sc_mi_array = np.loadtxt(nf1_sc_mi,delimiter=",")
nf2_sc_mi_array = np.loadtxt(nf2_sc_mi,delimiter=",")
nf4_sc_mi_array = np.loadtxt(nf4_sc_mi,delimiter=",")
nf8_sc_mi_array = np.loadtxt(nf8_sc_mi,delimiter=",")
nf16_sc_mi_array = np.loadtxt(nf16_sc_mi,delimiter=",")
nf32_sc_mi_array = np.loadtxt(nf32_sc_mi,delimiter=",")
# _____ Loading files End _____

# _____ Figure generation BEGIN _____
# designates colors to use
mdl_color = user_arguments["mdl_color"]
expt_color = user_arguments["expt_color"]

# obtains Cee-MI for different noise factors
nf1_Cee_MI = np.average(nf1_sc_mi_array)
nf2_Cee_MI = np.average(nf2_sc_mi_array)
nf4_Cee_MI = np.average(nf4_sc_mi_array)
nf8_Cee_MI = np.average(nf8_sc_mi_array)
nf16_Cee_MI = np.average(nf16_sc_mi_array)
nf32_Cee_MI = np.average(nf32_sc_mi_array)
nf_Cee_MI_list = [nf1_Cee_MI,nf2_Cee_MI,nf4_Cee_MI,nf8_Cee_MI,nf16_Cee_MI,nf32_Cee_MI]

# generates figure
fig,ax1 = plt.subplots(1)

# plots value of MI of CSAR
ax1.hlines(CSAR_mi_array,-1,6,linestyles="dashed",color = expt_color,linewidth = 3,label = "MI of Expt CSAR")

# plots Cee-MI values under different noise factors
ax1.bar(np.arange(6),nf_Cee_MI_list,color = mdl_color,label = "Model CeeMI")

# labels noise factors
ax1.set_xticks(np.arange(6))
ax1.set_xticklabels([1,2,4,8,16,32])

# sets limits
ax1.set_xlim([-.5,5.5])

# sets axis labels
ax1.set_xlabel("Noise Factor, $\eta$")
ax1.set_ylabel("Mutual Information (Bits)")

# removes spines
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

# adds legend
ax1.legend(ncol = 2,frameon=False, bbox_to_anchor=(.5,1.1),loc="upper center")

# generates and saves plot in the tight layout
plt.tight_layout()
plt.savefig(output_dir)
plt.savefig(f"{output_dir}.svg",format='svg')
plt.show()
# _____ Figure generation END _____
