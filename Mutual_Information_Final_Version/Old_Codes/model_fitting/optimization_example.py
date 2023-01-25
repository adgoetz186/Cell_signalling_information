from Mutual_Information_Final_Version.functions.Models.IGFR_Akt_FoxO import model
import numpy as np
from pysb.simulator import BngSimulator
import multiprocessing as mp
import functools
import scipy.optimize as so
import Mutual_Information_Final_Version.functions.file_processing_functions.user_input_functions as uif
import time
import random
import matplotlib.pyplot as plt
import copy as cp

plt.rcParams.update({'font.size': 15})


# This program generates outcomes of the IGFR/Akt/FoxO pathway using the Gillespie algorithm
# These results are used to test the moment closure accuracy


def sim_model(mdl_parameters, expt_values, times):
	akt_0 = mdl_parameters[12].astype(int)
	n_foxo_0 = (mdl_parameters[13] * mdl_parameters[10] / (mdl_parameters[10] + mdl_parameters[11])).astype(int)
	c_foxo_0 = (mdl_parameters[13] * mdl_parameters[11] / (mdl_parameters[10] + mdl_parameters[11])).astype(int) + 1
	igfr_0 = mdl_parameters[0].astype(int)
	new_params = np.zeros(np.size(mdl_parameters) + 3)
	new_params[0] = mdl_parameters[0] * mdl_parameters[1]
	new_params[1:12] = mdl_parameters[1:12]
	new_params[12] = igfr_0
	new_params[13] = akt_0
	new_params[14] = c_foxo_0
	new_params[15] = n_foxo_0
	error = 0
	for l in range(len(L)):
		new_params[16] = L[l]
		sim = BngSimulator(model, tspan=times, param_values=new_params)
		x = sim.run(n_runs=1, method='ode')
		dtf = x.dataframe
		error += np.sum((expt_values[l] - dtf["obs_FoxO_n"].to_numpy()) ** 2 / expt_values[l] ** 2)
	return error


# _____ Load arguments BEGIN _____
# Specifies path to file which specifies user specified arguments
user_input_file_path = "/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Figure_Generating_Programs/S_Gillespie_Supplimentary/General_Data_Generation/IGFR/ssa_IGFR_args.txt"
user_arguments = uif.read_file_to_arg_dict(user_input_file_path)
# _____ Load arguments END _____

# _____ File path declarations BEGIN _____
first_dose = "/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Data/IGFR/experimental_data/All_Doses_Properly_Weighted_With_Zero/0pm"
second_dose = "/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Data/IGFR/experimental_data/All_Doses_Properly_Weighted_With_Zero/10pm"
third_dose = "/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Data/IGFR/experimental_data/All_Doses_Properly_Weighted_With_Zero/15pm"
fourth_dose = "/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Data/IGFR/experimental_data/All_Doses_Properly_Weighted_With_Zero/20pm"
fifth_dose = "/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Data/IGFR/experimental_data/All_Doses_Properly_Weighted_With_Zero/25pm"
sixth_dose = "/Users/agoetz/PycharmProjects/Cell_signalling_information/Mutual_Information_Final_Version/Data/IGFR/experimental_data/All_Doses_Properly_Weighted_With_Zero/50pm"
# _____ File path declarations END _____

# _____ Loading files BEGIN _____
times = np.loadtxt(first_dose, delimiter=",")[0]
print(times)
first_dose = np.average(np.loadtxt(first_dose, delimiter=",")[1:], axis=0)
second_dose = np.average(np.loadtxt(second_dose, delimiter=",")[1:], axis=0)
fourth_dose = np.average(np.loadtxt(fourth_dose, delimiter=",")[1:], axis=0)
sixth_dose = np.average(np.loadtxt(sixth_dose, delimiter=",")[1:], axis=0)
titu = np.array([0, 1, 3, 5, 14, 19, 29, 59])
# fig,ax = plt.subplots(1)
# ax.scatter(times[titu],first_dose[titu],color = "black")
# ax.scatter(times[titu],second_dose[titu],color = "indigo")
# ax.scatter(times[titu],fourth_dose[titu],color = "green")
# ax.scatter(times[titu],sixth_dose[titu],color = "red")
# ax.plot(times[titu],first_dose[titu],color = "black")
# ax.plot(times[titu],second_dose[titu],color = "indigo")
# ax.plot(times[titu],fourth_dose[titu],color = "green")
# ax.plot(times[titu],sixth_dose[titu],color = "red")
# ax.spines["top"].set_visible(False)
# ax.spines["right"].set_visible(False)
# ax.set_ylabel("Nuclear Foxo")
# ax.set_xlabel("Time (Minutes)")
# plt.tight_layout()
# plt.show()
print(first_dose)
print(second_dose)
print(third_dose)
print(fourth_dose)
print(fifth_dose)
print(sixth_dose)

# _____ Loading files END _____

# _____ Main code BEGIN _____
# time allowed before sampling cells with no dose
t_f1 = 3600
t_f1_array = np.linspace(0, t_f1, 2)

# time points of main simulation
t_main = np.arange(0, 180 * 61, 180)
# list of doses
L = [0, 0.01, 0.02, 0.050]
dose_colors = ["black", "indigo", "green", "red"]
# number of cells to use
cell_count = user_arguments["cell_count"]

# The set of parameters used elsewhere doesnt allow for non equilibrium initial conditions
# this transforms the parameters into a set which can allow for a relevant subset of
# non equilibrium initial conditions
time_start = time.time()
# [0.9409369895991775, -2.7112006147516152, -1.124555325066475, -0.3412027350797875, -0.10028559912085329, 0.1287715298502503, -3.001369438724792, -5.753127968798565, -5.69776198643665, -3.338283637622721, -1.220670979074373, -2.7033120742251584, 3.652052848248105, 4.9870402869792665, 1.2304489213782739, 2.711807229041191, -inf]

# Bounds established from literature
mdl_lower_bounds = np.array([3, -3.5, -2, -4.5, -2, -1, -4, -6.5, -6.5, -4, -2, -3.5, 4.5, 2.8])
mdl_upper_bounds = np.array([4.5, -2, 0, -2, 0, 1, -2, -4.5, -4.5, -2, 0, -1.5, 5, 2.8])
mdl_middle_point = (mdl_lower_bounds + mdl_upper_bounds) / 2
mdl_point = 10 ** (
		mdl_lower_bounds + np.random.random(np.shape(mdl_upper_bounds)) * (mdl_upper_bounds - mdl_lower_bounds))


expt_doses = [first_dose[titu], second_dose[titu], fourth_dose[titu], sixth_dose[titu]]
times_to_sim = times[titu] * 60
bounded = so.Bounds(lb=10 ** mdl_lower_bounds, ub=10 ** mdl_upper_bounds)
res = so.minimize(sim_model, mdl_point, args=(expt_doses, times_to_sim), bounds=bounded, method='Nelder-Mead')
print(res)
akt_0 = res["x"][12].astype(int)
n_foxo_0 = (res["x"][13] * res["x"][10] / (res["x"][10] + res["x"][11])).astype(int)
c_foxo_0 = (res["x"][13] * res["x"][11] / (res["x"][10] + res["x"][11])).astype(int) + 1
igfr_0 = res["x"][0].astype(int)
new_params = np.zeros(np.size(res["x"]) + 3)
new_params[0] = res["x"][0] * res["x"][1]
new_params[1:12] = res["x"][1:12]
new_params[12] = igfr_0
new_params[13] = akt_0
new_params[14] = c_foxo_0
new_params[15] = n_foxo_0
colors = ["black", "indigo", "green", 'red']
fig, ax = plt.subplots(1)
for l in range(len(L)):
	new_params[16] = L[l]
	sim = BngSimulator(model, tspan=times_to_sim, param_values=new_params)
	x = sim.run(n_runs=1, method='ode')
	dtf = x.dataframe
	ax.plot(times_to_sim / 60, dtf["obs_FoxO_n"].to_numpy(), color=colors[l])
	ax.scatter(times_to_sim / 60, expt_doses[l], color=colors[l])
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_ylabel("Nuclear Foxo")
ax.set_xlabel("Time (Minutes)")
plt.tight_layout()
plt.show()