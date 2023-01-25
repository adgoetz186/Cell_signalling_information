import gillespy2
from gillespy2 import Model, Species, Parameter, Reaction
import numpy as np
import matplotlib.pyplot as plt
import testwhat
randlist = []
for i in range(8):
    randlist.append(testwhat.random())
input(randlist)

def runliststoaveragelist(rlists):
    averagelist = [0 for i in range(len(rlists[0]))]
    for i in rlists:
        for j in range(len(i)):
            averagelist[j] += i[j]
    for i in range(len(averagelist)):
        averagelist[i] /= len(rlists)
    print(averagelist)
    return averagelist

def stringify(timelist,valuelist):
    stringname = "{"
    for i in range(len(valuelist)):
        stringname += ("{" + str(timelist[i])+"," + str(valuelist[i]) + "},")
    return ((stringname + "}").replace("},}","}}"))

class Surfacereceptor(Model):
  def __init__(self, parameter_values = None):
    # Initialize the model.
    Model.__init__(self, name = "Surfacereceptor")

    # Define parameters.
    k1 = Parameter(name = 'k1', expression = 5)
    k2 = Parameter(name = 'k2', expression = 13)
    self.add_parameter([k1,k2])

    # Define molecular species.
    A = Species(name = 'A', initial_value = 10)
    self.add_species([A])

    # Define reactions.
    ab = Reaction(name = "ab", reactants = {A:1}, products = {},
                  propensity_function = "k1*A")
    bf = Reaction(name = "bf", reactants = {}, products = {A:1},
                  propensity_function = "k2")
    self.add_reaction([ab,bf])
    self.timespan(np.linspace(0, 2, 21))
model = Surfacereceptor()
ntraj = 10000
s_results = model.run(solver = gillespy2.solvers.NumPySSASolver,number_of_trajectories=ntraj)
Runaverage = s_results[0]
A1 = Runaverage['A']
A2 = Runaverage['A']**2
for index in range(1, ntraj):
    trajectory = s_results[index]
    A1+=trajectory["A"]
    A2 += np.square(trajectory["A"])
    #plt.plot(trajectory['time'], trajectory['A']**2, 'r')
    #plt.plot(trajectory['time'], trajectory['B']**2, 'b')
    #plt.plot(trajectory['time'], trajectory['F']**2, 'g')
    #plt.plot(trajectory['time'], trajectory['A']*trajectory['B'], 'r')
    #plt.plot(trajectory['time'], trajectory['A']*trajectory['F'], 'b')
    #plt.plot(trajectory['time'], trajectory['B']*trajectory['F'], 'g')
    #plt.show()
#plt.plot(Runaverage['time'], A2, 'r')
#plt.plot(Runaverage['time'], B2, 'b')
#plt.plot(Runaverage['time'], F2, 'g')
#plt.plot(Runaverage['time'], AB, 'r')
#plt.plot(Runaverage['time'], AF, 'b')
#plt.plot(Runaverage['time'], BF, 'g')
#plt.show()
A1/=ntraj
A2/=ntraj
print(stringify(Runaverage['time'],A1))
print(stringify(Runaverage['time'],A2))
plt.plot(Runaverage['time'],A1,'g')
plt.show()
plt.plot(Runaverage['time'],A2,'b')
plt.show()
input()
A1 = Runaverage['A'][-1]
B1 = Runaverage['B'][-1]
F1 = Runaverage['F'][-1]
AB = Runaverage['A'][-1]*Runaverage['B'][-1]
AF = Runaverage['A'][-1]*Runaverage['F'][-1]
BF = Runaverage['B'][-1]*Runaverage['F'][-1]
listAB = [AB-A1*B1]
listAF = [AF-A1*F1]
listBF = [BF-B1*F1]
listcount = [1]
for index in range(1, ntraj):
    trajectory = s_results[index]
    listcount.append(listcount[-1]+1)
    A1 += trajectory['A'][-1]
    B1 += trajectory['B'][-1]
    F1 += trajectory['F'][-1]
    AB += trajectory['A'][-1] * trajectory['B'][-1]
    AF += trajectory['A'][-1] * trajectory['F'][-1]
    BF += trajectory['B'][-1] * trajectory['F'][-1]
    listAB.append((AB/(index+1) - A1 * B1/(index+1)**2))
    listAF.append((AF/(index+1) - A1 * F1/(index+1)**2))
    listBF.append((BF/(index+1) - B1 * F1/(index+1)**2))
plt.plot(listcount,listAB,'r')
plt.plot(listcount,listAF,'b')
plt.plot(listcount,listBF,'g')
plt.show()
input([listAB[-1],listAF[-1],listBF[-1]])
covABlim = []
x = [i+1 for i in range(0,100)]
for index in range(0, 100):
    trajectory = s_results[index]
    varA += np.square(trajectory["A"] - averageA)
    covAB += np.multiply((trajectory["A"] - averageA),(trajectory["B"] - averageB))
    covABlim.append(covAB[-1]/(index+1))
    varB += np.square(trajectory["B"] - averageB)
    varF += np.square(trajectory["F"] - averageF)
plt.plot(x,covABlim)
plt.show()
varA/= 100
varB/= 100
varF/= 100
covAB/=100
plt.plot(Runaverage['time'], covAB, 'r')
plt.show()
plt.plot(Runaverage['time'], varA, 'r')
plt.plot(Runaverage['time'], varB, 'b')
plt.plot(Runaverage['time'], varF, 'g')
plt.show()
plt.plot(Runaverage['time'], averageA, 'r')
plt.plot(Runaverage['time'], averageB, 'b')
plt.plot(Runaverage['time'], averageF, 'g')
plt.show()
print(s_results)