from sympy import *
import numpy as np
from sympy.solvers.ode.systems import dsolve_system
t = Symbol("t")
#L,k1,k2,kn1,kn2,kpa,kda,kdf = symbols("L,k1,k2,kn1,kn2,kpa,kda,kdf")
L = 1
k1 = 2
k2 = 3
kn1 = 4
kn2 = 5
kpa = 6
kda = 7
kdf = 8
A = Function("A")(t)
B = Function("B")(t)
F = Function("F")(t)
A2 = Function("A2")(t)
B2 = Function("B2")(t)
F2 = Function("F2")(t)
AB = Function("AB")(t)
AF = Function("AF")(t)
BF = Function("BF")(t)
eqB2 =  2.0*AB*L*k1 + 1.0*A*L*k1 - 2.0*B2*k2 - 2.0*B2*kda - 2.0*B2*kn1 + 2.0*BF*kn2 + 1.0*B*k2 + 1.0*B*kda + 1.0*B*kn1 + 1.0*F*kn2 - diff(B2,t)
eqB = 1.0*A*L*k1 - 1.0*B*k2 - 1.0*B*kda - 1.0*B*kn1 + 1.0*F*kn2 -diff(B,t)
eqF = 1.0*B*k2 - 1.0*F*kdf - 1.0*F*kn2 - diff(F,t)
eqA = -1.0*A*L*k1 - 1.0*A*kda + 1.0*B*kn1 + 1.0*kpa -diff(A,t)
eqBF = 1.0*AF*L*k1 + 1.0*B2*k2 - 1.0*BF*k2 - 1.0*BF*kda - 1.0*BF*kdf - 1.0*BF*kn1 - 1.0*BF*kn2 - 1.0*B*k2 + 1.0*F2*kn2 - 1.0*F*kn2 -diff(BF,t)
eqAB = 1.0*A2*L*k1 - 1.0*AB*L*k1 - 1.0*AB*k2 - 2.0*AB*kda - 1.0*AB*kn1 + 1.0*AF*kn2 - 1.0*A*L*k1 + 1.0*B2*kn1 - 1.0*B*kn1 + 1.0*B*kpa - diff(AB,t)
eqF2 = 2.0*BF*k2 + 1.0*B*k2 - 2.0*F2*kdf - 2.0*F2*kn2 + 1.0*F*kdf + 1.0*F*kn2 - diff(F2,t)
eqAF = 1.0*AB*k2 - 1.0*AF*L*k1 - 1.0*AF*kda - 1.0*AF*kdf - 1.0*AF*kn2 + 1.0*BF*kn1 + 1.0*F*kpa -diff(AF)
eqA2 = -2.0*A2*L*k1 - 2.0*A2*kda + 2.0*AB*kn1 + 1.0*A*L*k1 + 1.0*A*kda + 2.0*A*kpa + 1.0*B*kn1 + 1.0*kpa -diff(A2)
print(dsolve_system([eqB,eqF,eqA],[A,B,F],t))
print(dsolve_system([eqB2,eqB,eqF,eqA,eqBF,eqAB,eqF2,eqAF,eqA2],[A,B,F,A2,B2,F2,AB,AF,BF],t))