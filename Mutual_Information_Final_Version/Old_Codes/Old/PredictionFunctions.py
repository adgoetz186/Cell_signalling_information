import numpy as np 
from scipy.integrate import solve_ivp
from MomentEquations import MomentsDiff_Eq_fn
from matplotlib import pyplot as plt


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
    #sol_dyn = solve_ivp(MomentsDiff_Eq_fn, tspan, z0, method='BDF', args=(K, IGF))
    sol_dyn = solve_ivp(MomentsDiff_Eq_fn, tspan, z0, method = 'RK23', args=(K,IGF))
    sol = sol_dyn.y[:,-1]
    return sol_dyn, sol

# Def the prediction funtion 
def FoxOn_preds_fn(k,L,t):
    """"Function spits out mean and second moment of nuclear foxO at time t
    input:  L = IGF concentration in nM 
            t = end time in seconds
    output: meanfoxO, variancefoxO"""
    _, Sol_t = solve_Moments_fn(k, L, t)
    foxOn_mean, foxOn_var = (Sol_t[7]), Sol_t[-1]
    # get the second moment from the variance 
    foxOn_s = foxOn_var + foxOn_mean**2
    return foxOn_mean, foxOn_s
    
    
def Moments_Preds_full_fn(k):
    """Function that gets the predictions for 6 concentrations each at 6 time points """
    nCons = 72 #number of constraints half of them means and half of them second moments 
    L  = np.array([10,15,20,25,50,25000])* 10 ** -3 #make it in nM
    times_arr = np.array([0,6,15,21,51,126])*60
    Moments_Preds_arr = np.zeros(nCons)
    i=0
    for igf in L: 
        for t in times_arr:
            Moments_Preds_arr[i], Moments_Preds_arr[i + int(nCons/2)] = FoxOn_preds_fn(k, igf, t)
            i+=1
    return Moments_Preds_arr
mod = 6
K = np.array([  -0.2, -4.1, -1.5, -4.8, -.5, -1.5,-2, -0.25 - mod,-5, -2.3, -2.7, -3.5, 4.2, 2.4])
K2 = np.array([ 0.2, -3.1, -0.5, -2.8,  .5, -.5, -.5, 1.25 - mod,-4, -1.8 ,-2.3, -2.8, 4.8, 3.1])
Kalt = np.array([ 0.2, -3.1, -.5,0.8010299956639813+1,  .5, -.5, -.5, 1.25 - 6,-5, -1.8 ,-2.3, -2.8, 4.8, 2.9])

data_matrix = np.loadtxt("params_22 (1).csv",delimiter=',')
print(data_matrix.shape)
input("d")
print(np.max(data_matrix[:,2]))
print(np.min(data_matrix[:,2]))
k_idx = np.random.randint(0,data_matrix.shape[0], 1)

data_matrix[:,2] -= 0.2


def display(moments_Preds_arr):
    L = np.array([10, 15, 20, 25, 50, 250]) * 10 ** -3  # make it in nM
    times_arr = np.array([0, 6, 15, 21, 51, 126]) * 60
    for i in range(len(moments_Preds_arr)):
        if i <=5:
            plt.plot(times_arr/3600,moments_Preds_arr[i],label=f"L = {L[i]*1000} pm")
    plt.legend()
    plt.ylabel("FoxOn")
    plt.xlabel("Time (hours)")
    plt.show()
    print(moments_Preds_arr)

#print(data_matrix[4575])
#input()
#display(np.reshape(Moments_Preds_full_fn(np.array([ 0.19965062, -3.17393893, -1.1812756,  -2.98212155,  0.22440149, -0.54112188,-0.54304167, -6.16375488, -4.11834815, -2.08577159, -2.64523385, -2.94373662,4.58544389,  2.4919829 ])),(12,-1)))
#input("pls stop")
#display(np.reshape(Moments_Preds_full_fn(Kalt),(12,-1)))
#input("test")

init = np.zeros_like(np.reshape(Moments_Preds_full_fn(data_matrix[k_idx[0]]),(12,-1)))
count = 0
for i in [175, 231, 484, 515, 893, 930, 1013, 1641, 1785, 2147, 2196, 2256, 2329, 2356, 2365, 2410, 2426, 2479, 2660, 2678, 2721, 2833, 2879, 2976, 3045, 3238, 3345, 3720, 4108, 4545, 4557, 4625, 5090, 5307]:
    init += np.reshape(Moments_Preds_full_fn(data_matrix[i]), (12, -1))
    count+=1
    print(count)

init/=count
display(init)

for i in range(len(data_matrix)):
    print(data_matrix[i].tolist())
    display(np.reshape(Moments_Preds_full_fn(data_matrix[i]), (12, -1)))

number_to_draw = 1
init = np.zeros_like(np.reshape(Moments_Preds_full_fn(data_matrix[k_idx[0]]),(12,-1)))
for i in k_idx:
    print(i)
    init+=np.reshape(Moments_Preds_full_fn(data_matrix[i]),(12,-1))
init/=number_to_draw
display(init)


