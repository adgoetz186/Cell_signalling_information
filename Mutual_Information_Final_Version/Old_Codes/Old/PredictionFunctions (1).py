import numpy as np 
from scipy.integrate import solve_ivp
from MomentEquations import MomentsDiff_Eq_fn
from matplotlib import pyplot as plt
import time

def Get_Init_fn(k): 
    """Inputs: k
    Output: z0 """
    # last two entries 
    nCom = 44
    k = 10**(k)

    k1 = k[0]   #Synthesis of IGFR 
    k2 = k[1]   #Degredation of IGFR
    k3 = k[2]   #Binding of IGFR to IGF
    k4 = k[3]   #Unbinding IGFR to IGF 
    k5 = k[4]   #Phosphorylation of bound receptor 
    k6 = k[5]   #Dephosphorylation of bound receptor 
    k7 = k[6]   #Dephosphorylation of AKT
    k8 = k[7]   #Phosphorylation of AKT
    k9 = k[8]   #Phosphorylation of FoxO
    k10 = k[9]  #Dephosphorylation of FoxO
    k11 = k[10] #Influx of FoxO to nucleus 
    k12 = k[11] #Efflux of FoxO from nucleus
    k_tot_Akt = k[12] 
    k_tot_foxo = k[13]
    # z[0]   #R
    # z[1]   #B
    # z[2]   #P
    # z[3]   #akt
    # z[4]   #pakt
    # z[5]   #pfoxoc
    # z[6]   #foxoc
    # z[7]   #foxon

    #
    z0 = np.zeros(nCom)
    # z0_1
    z0[0]= k1/k2
    # z0_3
    z0[3] = k_tot_Akt # 34050   # total AKT which can vary 
    # z0_7 
    z0[6] = (k_tot_foxo*k12)/(k11 + k12)   #710 is total FoxO which can also vary 
    # z0_8 
    z0[7] = (k_tot_foxo*k11)/(k11 + k12)  
    # z0_1_1
    z0[8]= k1/k2
    # z0_7_7
    # -(- 710*k11^2*k12 + 710*k11*k12^2)/(k11^2*(k11 + k12))
    z0[41] = (k_tot_foxo*k11*k12)/(k11**2 + 2*k11*k12 + k12**2)
    # z0_7_8 
    z0[42] =  -(k_tot_foxo*k11*k12)/((k11 + k12)**2)
    # z0_8_8
    z0[43] = (k_tot_foxo*k11*k12)/(k11**2 + 2*k11*k12 + k12**2)
    return z0

def solve_Moments_fn(K,IGF,tend):
    """Inputs: K , L = IGF , tend
    Outputs: Z solution"""
    tspan = np.array([0,tend])
    z0 = Get_Init_fn(K)
    sol_dyn = solve_ivp(MomentsDiff_Eq_fn, tspan, z0, method = 'LSODA', args=(K,IGF))
    #sol_dyn = solve_ivp(MomentsDiff_Eq_fn, tspan, z0, method='LSODA', args=(K, IGF))
    sol = sol_dyn.y[:,-1]
    return sol_dyn, sol

# Def the prediction funtion 
def FoxOn_preds_fn(k,L,t):
    """"Function spits out mean and second moment of nuclear foxO at time t
    input:  L = IGF concentration in nM 
            t = end time in seconds
    output: meanfoxO, variancefoxO"""
    _, Sol_t = solve_Moments_fn(k, L, t)
    foxOn_mean, foxOn_var = Sol_t[7], Sol_t[-1]
    # get the second moment from the variance 
    foxOn_s = foxOn_var + foxOn_mean**2
    return foxOn_mean, foxOn_s
    
    
def Moments_Preds_full_fn(k):
    """Function that gets the predictions for 6 concentrations each at 6 time points """
    nCons = 72 #number of constraints half of them means and half of them second moments 
    L  = np.array([10,15,20,25,50,250]) * 10 ** -3 #make it in nM
    times_arr = np.array([0,6,15,21,51,126])*60
    Moments_Preds_arr = np.zeros(nCons)
    i=0
    for igf in L: 
        for t in times_arr:
            Moments_Preds_arr[i], Moments_Preds_arr[i + int(nCons/2)] = FoxOn_preds_fn(k, igf, t)
            i+=1

    return Moments_Preds_arr

def display(moments_Preds_arr):
    L = np.array([10, 15, 20, 25, 50, 250]) * 10 ** -3  # make it in nM
    times_arr = np.array([0, 6, 15, 21, 51, 126]) * 60
    for i in range(len(moments_Preds_arr)):
        if i <=5:
            plt.plot(times_arr/3600,moments_Preds_arr[i],label=f"L = {L[i]}")
    plt.legend()
    plt.ylabel("FoxOn")
    plt.xlabel("Time (hours)")
    plt.title("Scipy BDF solver")
    plt.show()


data_matrix = np.loadtxt("params_22 (1).csv",delimiter=',')
data_matrix[:,2] -= 0.2
time_start = time.time()
error_index = []
for i in range(len(data_matrix)):
    print(i)
    test_feature = np.reshape(Moments_Preds_full_fn(data_matrix[i]), (12, -1))
    eq = True
    for j in range(1,6):
        if np.any(test_feature[j]>test_feature[j-1]):
            eq = False
    if eq == False:
        error_index.append(i)
        print(i)
        print(data_matrix[i].tolist())
        print("lookup")
    if i%100 == 0:
        print(f"Took {time.time()-time_start} seconds to do {i} runs")
print(error_index)
print(f"Took {time.time()-time_start} seconds to do {len(data_matrix)} runs")
#display(np.reshape(Moments_Preds_full_fn(data_matrix[i]), (12, -1)))
