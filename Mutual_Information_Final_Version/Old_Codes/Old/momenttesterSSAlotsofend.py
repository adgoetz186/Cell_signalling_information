import gillespy2
from gillespy2 import Model, Species, Parameter, Reaction
import numpy as np
import matplotlib.pyplot as plt

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
    k1 = Parameter(name = 'k1', expression = 1)
    k2 = Parameter(name = 'k2', expression = 2)
    kn1   = Parameter(name = 'kn1',   expression = 3)
    kn2  = Parameter(name = 'kn2',  expression = 4)
    kpa     = Parameter(name = 'kpa',     expression = 5)
    kda = Parameter(name='kda', expression=6)
    kdf = Parameter(name='kdf', expression=7)
    L = Parameter(name='L', expression=8)
    self.add_parameter([k1,k2,kn1,kn2,kpa,kda,kdf,L])

    # Define molecular species.
    A = Species(name = 'A', initial_value = 0)
    B = Species(name = 'B', initial_value = 0)
    F = Species(name='F', initial_value=0)
    self.add_species([A,B,F])

    # Define reactions.
    ab = Reaction(name = "ab", reactants = {A:1}, products = {B:1},
                  propensity_function = "k1*A*L")
    bf = Reaction(name = "bf", reactants = {B:1}, products = {F:1},
                  propensity_function = "k2*B")
    ba = Reaction(name = "ba", reactants = {B:1}, products = {A:1},
                  propensity_function = "kn1*B")
    fb = Reaction(name = "fb", reactants = {F:1}, products = {B:1},
                  propensity_function = "kn2*F")

    pa = Reaction(name="pa", reactants={}, products={A: 1},
                  propensity_function="kpa")
    da = Reaction(name="da", reactants={A: 1}, products={},
                  propensity_function="kda*A")
    db = Reaction(name="db", reactants={B: 1}, products={},
                  propensity_function="kda*B")
    df = Reaction(name="df", reactants={F: 1}, products={},
                  propensity_function="kdf*F")
    self.add_reaction([ab,bf,fb,ba,pa,da,db,df])
    self.timespan(np.linspace(0, 2, 21))
bleh = []
minil = []
avgl = []
for i in range(1000):
    model = Surfacereceptor()
    ntraj = 1000
    s_results = model.run(solver = gillespy2.solvers.NumPySSASolver,number_of_trajectories=ntraj)
    Runaverage = s_results[0]
    A1 = Runaverage['A']
    A2 = Runaverage['A']**2
    B1 = Runaverage['B']
    B2 = Runaverage['B']**2
    F1 = Runaverage['F']
    F2 = Runaverage['F']**2
    AB = Runaverage['A']*Runaverage['B']
    AF = Runaverage['A']*Runaverage['F']
    BF = Runaverage['B']*Runaverage['F']
    for index in range(1, ntraj):
        trajectory = s_results[index]
        A1+=trajectory["A"]
        A2 += np.square(trajectory["A"])
        B1 += trajectory["B"]
        B2 += np.square(trajectory["B"])
        F1 += trajectory["F"]
        F2 += np.square(trajectory["F"])
        AB += np.multiply(trajectory['A'],trajectory['B'])
        AF += np.multiply(trajectory['A'],trajectory['F'])
        BF += np.multiply(trajectory['B'] , trajectory['F'])
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
    B1/=ntraj
    F1/=ntraj
    A2/=ntraj
    B2/=ntraj
    F2/=ntraj
    AB/=ntraj
    AF/=ntraj
    BF/=ntraj

    R1 = A1+B1+F1
    R2 = A2+B2+F2+2*AB+2*AF+2*BF
    covAB = AB-A1*B1
    covAF = AF-A1*F1
    covBF = BF-B1*F1
    print(stringify(Runaverage['time'],R1))
    print(stringify(Runaverage['time'],R2))
    print(stringify(Runaverage['time'],covAB))
    print(stringify(Runaverage['time'],covAF))
    print(stringify(Runaverage['time'],covBF))
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
    #plt.show()
    print([listAB[-1],listAF[-1],listBF[-1]])
    bleh.append([listAB[-1],listAF[-1],listBF[-1]])
    mini = [1,1,1]
    avg = [0,0,0]
    for i in bleh:
        for j in range(3):
            mini[j] = min(mini[j],abs(i[j]))
            avg[j]+=i[j]
    for i in range(len(avg)):
        avg[i] /= len(bleh)
    minil.append(mini)
    avgl.append(avg)
plt.show()
num = [i for i in range(1000)]
plt.plot(num,minil,'b')
plt.plot(num,avgl,'g')
plt.show()