import numpy as np
import sys
import os
from pathlib import Path, PureWindowsPath
import Mutual_Information_Main.functions.data_processing_functions.params_to_moments_igf as p2m

# Generates moments of responses for IGF/Akt/FOXO pathway for comparison with the step dose experimental data

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

# _____ File path declarations BEGIN _____
# Specifies path to the parameter draws obtained from MERIDIAN inference
params_path = Path("Data/IGFR/Model_Params/params_162/params_162.csv")

# Specifies path to output
output_path = Path("Data/IGFR/Moments/Model_Moments/4_dose_response/moments_162_4_dose_60_90_min")

# Specifies path to output header
output_path_header = Path("Data/IGFR/Moments/Model_Moments/4_dose_response/moments_162_4_dose_60_90_min_header_structure.csv")
# _____ File path declarations END _____

# _____ Loading files BEGIN _____
# Loads parameter draws from MERIDIAN inference
params = np.loadtxt(params_path,delimiter=",")
# _____ Loading files END _____


# _____ Main code BEGIN _____
# number of cells to generate (50000)
number_of_cells_to_use = 50000
L = [0,17.5,37.5,125]
L_names = ["0_pm","17.5_pm","37.5_pm","125_pm"]
times_arr = np.array([60, 90])*60
times_names = [f"{int(i/60)}_min" for i in times_arr]
L = [i/1000 for i in L]
nd = len(L)
nt = len(times_names)

# creates list of cell numbers which will be used
cell_numbers_to_use = np.arange(np.shape(params)[0])[:number_of_cells_to_use]

# Creates array of responses to given doses
# Ordered in the following way
# Moment 1 Dose 1 Time 1, ..., Moment 1 Dose 1 Time n, ..., Moment 1 Dose n Time n, ..., Moment 2 Dose n Time n
dose_response_array = np.zeros((number_of_cells_to_use,int(nd*2*nt)))
count = 0
for i in cell_numbers_to_use:
	response = p2m.cell_pred_fn(params[i], times_arr, L, meth='BDF')
	means = response[0]
	var = response[1]
	scnd_mom = means**2+var
	dose_response_array[count,:int(nd*nt)] = means
	dose_response_array[count, int(nd*nt):] = scnd_mom
	if count % 10 == 0:
		print((count+1)/number_of_cells_to_use)
	count+=1
np.save(output_path,dose_response_array)

# Creates headers for each column of the response array
# Ordered in the following way
# Moment 1 Dose 1 Time 1, ..., Moment 1 Dose 1 Time n, ..., Moment 1 Dose n Time n, ..., Moment 2 Dose n Time n
moment_list = ["first","second"]
header_list = []
for moment in moment_list:
	for dose in L_names:
		for time in times_names:
			header_list.append({"moment":moment,"dose":dose,"time":time})
with open(output_path_header,'w') as output_header:
	for i in header_list:
		output_header.write(str(i)+"\n")
# _____ Main code END _____
