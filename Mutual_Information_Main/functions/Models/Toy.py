from pysb import *

# Parameters in order:
# kprod - unbound unphosphorylated receptor synthesis rate
# kbind - receptor ligand binding rate constant
# kunbind - receptor ligand unbinding rate constant
# kdeg - non activated receptor degredation rate
# R_0 - Initial unbound R count
# L - Ligand concentration


Model()
Compartment('e',dimension=2,parent=None)

# Due to the replacement rate of the ligand, the ligand concentration acts like a constant
# so IGF binding acts as a change of state of IGFR rather than a binding event
Monomer('R', ['ec'],{'ec':['b','ub']})

Parameter('kprod',1.0)
Parameter('kdeg',1.0)
Parameter('kbind',1.0)
Parameter('kunbind',1.0)
Parameter('R_0', 1)
Parameter('L', 1.0)

Expression("keffbind",kbind*L)

Initial(R(ec='ub')**e, R_0)

Rule('r_prod', None >> R(ec='ub')**e,kprod)
Rule('r_deg', R()**e >> None ,kdeg)
Rule('r_bind', R(ec='ub')**e | R(ec='b')**e,keffbind,kunbind)



Observable('obs_R', R(ec='ub'))
Observable('obs_B', R(ec='b'))
Observable('obs_R_totalB', R())
