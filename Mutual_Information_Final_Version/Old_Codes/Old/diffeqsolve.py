import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import math

# function that returns dy/dt
def model(U,t):
    Nmkr, mRNA = U
    kmprod = math.log(2) / ((8.25) * (60 ** 2)) * 10.44
    kmdeg = math.log(2) / ((8.25) * (60 ** 2))
    knprod = 21.58 / (60 ** 2)
    kndeg = math.log(2) / ((22.23) * (60 ** 2))
    dmRNAdt = kmprod-kmdeg*mRNA
    dnmkrdt = knprod*Nmkr-kndeg*Nmkr
    return [dmRNAdt,dnmkrdt]

# initial condition
y0 = 5

# time points
t = np.linspace(0,400000,200)

# solve ODE
y = odeint(model,[0,0],t)

# plot results
plt.plot(t,y)
plt.xlabel('time')
plt.ylabel('y(t)')
plt.show()