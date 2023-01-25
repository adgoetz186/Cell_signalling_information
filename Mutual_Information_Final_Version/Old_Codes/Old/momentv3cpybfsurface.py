from sympy import *
import numpy as np




x = Symbol('x')
f = Function('f')

print(f)
print(f(3*x))
names = ["A","B","C","D"]

#Enter explicit values here:
#reactionvalues = {"k1":1,"k2":2,"kn1":3,"kn2":4,"kpa":5,"kda":6,"kdf":7,"L":8}
reactionvalues = {"k3":0.01,"k5":0.011,"k4":0.005,"k6":0.00022,"V":1.66 * 10 ** (-13)}

k1 = ["A->B","A*k1*tf*V"]
k2 = ["B->A","B*k2"]
k3 = ["0->C","B*k3"]
k4 = ["C->0","C*k4"]
k5 = ["0->D","C*k5"]
k6 = ["D->0","D*k6"]


rlist = [k1,k2,k3,k4,k5,k6]
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


todolist = ["D**2"]
done = []
eqlist = []
constants = []
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
D1,D2 = symbols("D1,D2")
for keys in solution[0].keys():
    print(keys, solution[0][keys])
print(solution[0][D1])
print(solution[0][D2]-solution[0][D1]**2)
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


