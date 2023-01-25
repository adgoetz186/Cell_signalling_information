import numpy as np
import Mutual_Information_Final_Version.functions.file_processing_functions.user_input_functions as uif
import Mutual_Information_Final_Version.functions.data_processing_functions.params_to_moments_igf as p2m

# Generates moments of responses for IGF/Akt/FOXO pathway

# _____ Load arguments BEGIN _____
# Specifies path to file which specifies user specified arguments for mapping cell parameters to cell responses
user_input_file_path = "/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Figure_Generating_Programs/Figure_3/General_Data_Generation/igfr_params_to_moments/3_dose_args.txt"
user_arguments = uif.read_file_to_arg_dict(user_input_file_path)
# _____ Load arguments END _____


# _____ Loading files BEGIN _____
# Loads parameter draws from MERIDIAN inference
params = np.loadtxt(user_arguments["parameter_file_path"],delimiter=",")
# _____ Loading files END _____


# _____ Main code BEGIN _____
number_of_cells_to_use = user_arguments["number_of_cells_to_use"]
L = user_arguments["L"]
L_names = user_arguments["L_names"]
times_arr = np.array(user_arguments["times"])*60
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
np.savetxt(user_arguments["output_file_path"],dose_response_array,delimiter=",")

# Creates headers for each column of the response array
# Ordered in the following way
# Moment 1 Dose 1 Time 1, ..., Moment 1 Dose 1 Time n, ..., Moment 1 Dose n Time n, ..., Moment 2 Dose n Time n
moment_list = ["first","second"]
header_list = []
for moment in moment_list:
	for dose in L_names:
		for time in times_names:
			header_list.append({"moment":moment,"dose":dose,"time":time})
with open(user_arguments["output_header_file_path"],'w') as output_header:
	for i in header_list:
		output_header.write(str(i)+"\n")
# _____ Main code END _____
