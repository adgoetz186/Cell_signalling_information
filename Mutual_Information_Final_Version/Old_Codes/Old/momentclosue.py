from sympy import *
import numpy as np




x = Symbol('x')
f = Function('f')

print(f)
print(f(3*x))
names = ["A","B","F","G","H","K","T","V"]
reactionvalues = {"kn1":10**(-0.52),"k2":10**(-6.3),"kn2":10**(-1.25),"kp":10**0.05,"kdp":10**(-1.15),"kprod":16.539,"krec":10**(-3.55),"ksrec":10**(-4.42),"kdeg":10**(-3.6),"ksdeg":10**(-2.5)}

kap = ["0->A","kap"]

k1 = ["A->B","L*k1*A"]
k2 = ["B->F","k2*B"]
kn1 = ["B->A","kn1*B"]
kn2 = ["F->B","kn2*F"]
kp = ["F->G","kp*F"]
kdp = ["G->F","kdp*G"]

kai = ["A->H","A*ki"]
karec = ["H->A","H*krec"]

kbi = ["B->K","B*ki"]
kbrec = ["K->B","K*krec"]

kfi = ["F->T","F*ki"]
kfrec = ["T->F","T*krec"]

kgi = ["G->V","G*ksi"]
kgrec = ["V->G","V*ksrec"]

khdeg = ["H->0","kdeg*H"]
kkdeg = ["K->0","kdeg*K"]
ktdeg = ["T->0","kdeg*T"]
kvdeg = ["V->0","ksdeg*V"]

kpi = ["T->V","kp*T"]
kdpi = ["V->T","kdp*V"]

rlist = [kap,k1,k2,kn1,kn2,kp,kdp,kai,karec,kbi,kbrec,kfi,kfrec,kgi,kgrec,khdeg,kkdeg,ktdeg,kvdeg,kpi,kdpi]
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
    print(expv)
    input(listofvalues)

    for i in listofvalues:
        print(i.args)
        if len(i.args) <= 1:
            print([sympynames,i])
            u[sympynames.index(i)] = 1
        else:
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
        #value = value.subs(sympify("Off"), sympify("1-On"))
        #value = value.subs(sympify("L"), sympify("50-On"))
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
    for k in range(len(names)):
        if degree(polyentry, sympify(names[k])) > 0:
            namecomp += (names[k] + str(degree(polyentry, sympify(names[k]))))
    return namecomp


todolist = ["G**2"]
done = []
eqlist = []
while len(todolist) >= 1:
    expvalue = expectation(rlist,names, todolist[0],limit = 2)
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
    del todolist[0]
    print(eqentry[1])
    if eqentry[1] != 0:
        for i in eqentry[1].args:
            addtotodo = 1
            for j in range(1,len(i.args)):
                print(i.args[j])
                addtotodo *= i.args[j]
            if addtotodo != 1 and addtotodo not in done and addtotodo not in todolist:
                todolist.append(addtotodo)
        eqlist.append(eqentry)
print(len(eqlist))



eqtosolve = []
for i in range(len(eqlist)):
    firstentry = sympify(polyentrytosymstring(eqlist[i][0]))
    equation = eqlist[i][1]
    buildterm = 0
    while equation != 0:
        if equation.is_constant():
            buildterm += equation
            equation = 0
        else:
            ft = LT(equation)
            equation -= ft
            buildterm += (sympify(polyentrytosymstring(ft)) * LC(ft))
    eqtosolve.append([firstentry,buildterm])
solver = [eqtosolve[i][1] for i in range(len(eqtosolve))]
solution = solve(solver,set = True)
print(solution)



