from sympy import *
import numpy as np
from scipy.stats import nbinom
import seaborn as sns
from matplotlib import pyplot as plt
import time
import matplotlib.ticker as ticker


x = Symbol('x')
f = Function('f')

print(f)
print(f(3*x))
names = ["A","B"]

#Enter explicit values here:
#reactionvalues = {"k1":1,"k2":2,"kn1":3,"kn2":4,"kpa":5,"kda":6,"kdf":7,"L":8}
reactionvalues = {"k3":0.01,"k5":0.011,"k4":0.005,"k6":0.00022,"V":1.66 * 10 ** (-13)}

k1 = ["A->B","A*k1*L"]
k2 = ["B->A","B*k2"]



rlist = [k1,k2]
for key in reactionvalues.keys():
    for i in range(len(rlist)):
        rlist[i][1] = rlist[i][1].replace(key,str(reactionvalues[key]))
print(rlist)
def a(reactionlist,names):
    a = np.zeros([len(names), len(reactionlist)])
    for j in range(len(reactionlist)):
        print(reactionlist[j][0])
        print(reactionlist[j][0].split("->"))
        print(reactionlist[j][0].split("->")[0])
        print(reactionlist[j][0].split("->")[0].split("+"))
        reactants = reactionlist[j][0].split("->")[0].split("+")
        for i in range(len(reactants)):
            if "*" in reactants[i]:
                reactants[i] = [reactants[i].split("*")[1],reactants[i].split("*")[0]]
            elif reactants[i] == "0":
                reactants[i] = [names[0], "0"]
            else:
                reactants[i] = [reactants[i],"1"]
        for i in reactants:
            indexofreact = names.index(i[0])
            a[indexofreact,j] += eval(i[1])
    print(a)
    return a

def b(reactionlist,names):
    b = np.zeros([len(names), len(reactionlist)])
    for j in range(len(reactionlist)):
        products = reactionlist[j][0].split("->")[1].split("+")
        for i in range(len(products)):
            if "*" in products[i]:
                products[i] = [products[i].split("*")[1],products[i].split("*")[0]]
            elif products[i] == "0":
                products[i] = [names[0], "0"]
            else:
                products[i] = [products[i],"1"]
        for i in products:
            indexofprod = names.index(i[0])
            b[indexofprod,j] += eval(i[1])
    print(b)
    return b
Amat = a(rlist,names)
Bmat = b(rlist,names)


def gamma(amat,bmat):
    print(bmat-amat)
    return bmat-amat
Gmat = gamma(Amat,Bmat)
gammalist = []
for i in range(len(rlist)):
    gammalist.append(Gmat[:,i])
print(gammalist)

def M(names,u):
    value = 1
    for i in range(len(names)):
        evaluate = names[i] + "**" + str(u[i])
        symbolicexp = sympify(evaluate)
        print(symbolicexp)
        value *= symbolicexp
    print(expand(value))
    return value


def expressionTou(names,expv):
    if expv.func == Pow or expv.func == Symbol:
        listofvalues = (expv,)
    else:
        listofvalues = expv.args
    u = [0 for i in range(len(names))]
    sympynames = [sympify(names[i]) for i in range(len(names))]
    print(len(listofvalues))


    for i in listofvalues:
        print(i.args)
        if len(i.args) <= 1:
            if i in sympynames:
                print([sympynames,i])
                u[sympynames.index(i)] = 1
        else:
            if i.args[0] in sympynames:
                u[sympynames.index(i.args[0])] = i.args[1]
    print(u)

    return u

def expectation(reactionlist,names,expv,**kwargs):
    if type(expv) is str:
        expv = sympify(expv)
    u = expressionTou(names,expv)
    Amat = a(reactionlist, names)
    Bmat = b(reactionlist, names)
    G = gamma(Amat,Bmat)
    eq = 0
    for i in range(len(reactionlist)):
        name2 = ["(" + str(names[j]) + "+" + str(G[j,i]) + ")" for j in range(len(names))]
        print(name2)
        activation = sympify(reactionlist[i][1])
        value = (M(name2,u) - M(names,u))*activation
        #Adding assumptions here seems the best way mathmatically
        #Cosider f(x,y) = xy+x^2y^2 if we know y = 1-x
        #and want to only keep second order terms, we want to
        #substitute before we eliminate higher order terms
        #thats what this approach does
        #f(x,y) = x + O(3) under this assumption
        #value = value.subs(sympify("A"), sympify("r-B-F"))
        value = value.subs(sympify("A"), sympify("1-B"))
        if "limit" in kwargs.keys():
            valuebuild = 0
            if value != 0:
                valuebuffer = poly(expand(value))
                while valuebuffer != 0:
                    valuetoconsider = LT(valuebuffer)
                    valuebuffer -= valuetoconsider
                    order = 0
                    for j in names:
                        order += degree(valuetoconsider,sympify(j))
                    if order <= kwargs["limit"]:
                        valuebuild += valuetoconsider
                    else:
                        valuebuild += 0
            value = valuebuild
        #Make all basic terms into constants On*Off
        eq += expand(value)
    print(eq)
    if "limit" in kwargs.keys():
        if sum(u) > kwargs["limit"]:
            eq = 0
    return eq

def polyentrytosymstring(polyentry):
    namecomp = ""
    symbolvalue = 1
    for k in range(len(names)):
        if degree(polyentry, sympify(names[k])) > 0:
            namecomp += (names[k] + str(degree(polyentry, sympify(names[k]))))
            symbolvalue *= sympify(names[k] + "**" + str(degree(polyentry, sympify(names[k]))))
    return [namecomp,polyentry/symbolvalue]


todolist = ["A**2"]
done = []
eqlist = []
constants = ["A-B=0"]
sympynames = [sympify(names[i]) for i in range(len(names))]
while len(todolist) >= 1:
    expvalue = expectation(rlist,names, todolist[0])
    expvalue = expand(expvalue)
    # Adding assumptions here seems flawed, but might be the way to go
    # Check if other papers do this, if they do then they are likely wrong
    # Should be provable that other way converges to solution faster
    # Cosider f(x,y) = xy+x^2y^2 if we know y = 1-x
    # and want to only keep second order terms,
    # this approach is like cutting off all higher order terms before substituting
    # f(x,y) = x + x^2 + O(3) under this assumption
    #expvalue = expvalue.subs(sympify("Off"), sympify("1-On"))
    #expvalue = expand(expvalue.subs(sympify("L"), sympify("50-On")))
    eqentry = [sympify(todolist[0]),expvalue]
    done.append(sympify(todolist[0]))
    #input(todolist)
    del todolist[0]
    if eqentry[1] != 0:
        for i in eqentry[1].args:
            addtotodo = 1
            for j in range(0,len(i.args)):
                print(i.args[j])
                dertest = 0
                for k in sympynames:
                    dertest += (diff(i.args[j],k))**2
                if dertest != 0:
                    addtotodo *= i.args[j]
            if addtotodo != 1 and addtotodo not in done and addtotodo not in todolist:
                todolist.append(addtotodo)
        eqlist.append(eqentry)
print(len(eqlist))
for i in eqlist:
    print(i)
#input()
print("ok")


eqtosolve = []
for i in range(len(eqlist)):
    firstentry = sympify(polyentrytosymstring(eqlist[i][0])[0])
    equation = eqlist[i][1]
    buildterm = 0
    while equation != 0:
        if equation.is_constant():
            buildterm += equation
            equation = 0
        else:
            ft = LT(equation)
            equation -= ft
            symstringlist = polyentrytosymstring(ft)
            if symstringlist[0] != "":
                buildterm += (sympify(symstringlist[0])) * sympify(symstringlist[1])
            else:
                buildterm += sympify(symstringlist[1])
    eqtosolve.append([firstentry,buildterm])
print(eqtosolve)
solver = [eqtosolve[i][1] for i in range(len(eqtosolve))]
tosolve = [eqtosolve[i][0] for i in range(len(eqtosolve))]
print(solver)
print(tosolve)
solution = solve(solver,tosolve,dict = True)
print(solution)
print(solution[0])
input("solution")
D1,D2 = symbols("A1,A2")
for keys in solution[0].keys():
    print(keys, solution[0][keys])
print(solution[0][D1])
print(solution[0][D2]-solution[0][D1]**2)
input("test")
mean = solution[0][D1]
variance = solution[0][D2]-solution[0][D1]**2

scale = variance/mean
shape = mean/scale
print("mean")
print(mean)
print(variance)
print("ok")
print(scale)
print(shape)
response_min = 0
response_max = 200
input_list = [0,1,10,100,1000,10000]
number_of_inputs = len(input_list)
kstep = 0.05
response_values = np.arange(response_max-response_min+1)
print(len(response_values))
column_count_heatmap = 160
row_count_heatmap = 100
Kalist = [10**(3+i*kstep) for i in range(column_count_heatmap)]
Kblist = [10**(-5+i*kstep) for i in range(row_count_heatmap)]
mutual_info_matrix = np.empty((len(Kblist),len(Kalist)))
start_time = time.time()
for k1_index in range(len(Kalist)):
    print(k1_index)
    for k2_index in range(len(Kblist)):
        k1 = Kalist[k1_index]
        k2 = Kblist[k2_index]
        conditional_response = np.empty((number_of_inputs, response_max - response_min + 1))
        for tf_index in range(number_of_inputs):
            tf = input_list[tf_index]
            if tf > 0:
                meanval = N(mean.subs(sympify("k1"), k1).subs(sympify("k2"), k2).subs(sympify("tf"), tf))
                varianceval = N(variance.subs(sympify("k1"), k1).subs(sympify("k2"), k2).subs(sympify("tf"), tf))
                p = meanval / varianceval
                n = meanval * p / (1 - p)
                #input(nprime)
                #nbinom.pmf(response_values, n, 0.5)
                #nbinom.pmf(response_values, 10, p)
                #nbinom.pmf(1, n, p)
                #nbinom.pmf(response_values, float(n), float(p))

                conditional_response[tf_index] = nbinom.pmf(response_values,  float(n), float(p))
            else:
                conditional_response[tf_index] = np.append(1,np.zeros(response_max - response_min))
        def mutual_information_from_matrix(relative_input_vector):
            # print(f"Input vector is {relative_input_vector}")
            signal_ratios = [relative_input_vector]
            response_values = np.matmul(signal_ratios, conditional_response)
            signal_matrix = np.transpose(np.repeat(signal_ratios, repeats=np.shape(conditional_response)[1], axis=0))
            response_matrix = np.repeat(response_values, repeats=np.shape(conditional_response)[0], axis=0)
            ratio = np.divide(conditional_response, response_matrix, out=np.zeros_like(conditional_response),
                              where=response_matrix != 0)
            log_of_ratio = np.log2(ratio, out=np.zeros_like(ratio), where=ratio != 0)
            mutual_info_matrix = np.multiply(np.multiply(signal_matrix, conditional_response), log_of_ratio)
            return np.sum(np.nan_to_num(mutual_info_matrix, nan=0))

        uniform_input_vector = np.ones(number_of_inputs)/number_of_inputs
        mi = mutual_information_from_matrix(uniform_input_vector)
        mutual_info_matrix[k2_index,k1_index] = mi
        #ax = sns.heatmap(conditional_response)
        #plt.show()
print(time.time()-start_time)
#ax = sns.heatmap(mutual_info_matrix,cmap=sns.color_palette("viridis", as_cmap=True))
ax = sns.heatmap(mutual_info_matrix,cmap=sns.color_palette("Paired"))

ax.xaxis.set_major_locator(plt.MaxNLocator(4))
tickstart = np.round(np.log10(Kalist),decimals=1)[0:4]
print(tickstart)
values = np.arange(4)*3
print(values)
print(values+tickstart)
print(np.round(np.log10(Kalist),decimals=1))
all_tick = np.round(np.log10(Kalist))[0:len(Kalist):2]
print(all_tick)
#input()
tick_to_use = []
count = 0
for i in all_tick:
    if count%3==0:
        tick_to_use.append(i)
    count+=1
print(tick_to_use)

ticks_to_use = np.arange(0,column_count_heatmap,20)+.5
labels_to_use_base = np.arange(0,column_count_heatmap,20)
print(tick_to_use)
lables_to_use = (labels_to_use_base+6)/2
plt.xticks(ticks_to_use,labels=lables_to_use)


ticks_to_use = np.arange(0,row_count_heatmap,20)+.5
labels_to_use_base = np.arange(0,row_count_heatmap,20)
lables_to_use = (labels_to_use_base-10)/2
plt.yticks(ticks_to_use,labels=lables_to_use)


#ax.xaxis.set_major_locator(ticker.MultipleLocator(2))

plt.title("Information")
plt.xlabel("$log_{10}k_{on}s^{-1}M^{-1}$")
plt.ylabel("$log_{10}k_{off}s^{-1}$")
plt.show()

