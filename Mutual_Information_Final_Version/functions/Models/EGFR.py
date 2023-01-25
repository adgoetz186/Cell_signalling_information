from pysb import *

# Parameters in order:
# kprod - unbound unphosphorylated receptor synthesis rate
# kbind - receptor ligand binding rate constant
# kunbind - receptor ligand unbinding rate constant
# kp - receptor phosphorylation rate constant
# kdp - receptor dephosphorylation rate constant
# kdeg - non activated receptor degredation rate
# ksdeg - activated receptor degredation rate
# EGFR_0 - Initial unphosphorylated unbound EGFR count


Model()
Compartment('e',dimension=2,parent=None)

# Due to the replacement rate of the ligand, the ligand concentration acts like a constant
# so IGF binding acts as a change of state of IGFR rather than a binding event
Monomer('EGFR', ['k','ec'],{'k':['up','p'],'ec':['b','ub']})

Parameter('kprod',1.0)
#print(kprod.value)
# this is effective binding constant, should be igf concentration and rate constant
Parameter('kbind',1.0)
Parameter('kunbind',1.0)
Parameter('krp',1.0)
Parameter('krdp',1.0)
Parameter('kdeg',1.0)
Parameter('ksdeg',1.0)

Parameter('EGFR_0', 1)

Parameter('EGF', 1.0)

Expression("keffbind",kbind*EGF)

Initial(EGFR(k='up',ec='ub')**e, EGFR_0)

Rule('r_prod', None >> EGFR(k='up',ec='ub')**e,kprod)
Rule('r_deg', EGFR(k='up')**e >> None ,kdeg)
Rule('r_sdeg', EGFR(k='p')**e >> None ,ksdeg)
Rule('r_bind', EGFR(ec='ub',k='up')**e | EGFR(ec='b',k='up')**e,keffbind,kunbind)
Rule('rb_phos', EGFR(k='up',ec='b')**e | EGFR(k='p',ec='b')**e,krp,krdp)



Observable('obs_EGFR', EGFR(k='up',ec='ub'))
Observable('obs_EGFRb', EGFR(k='up',ec='b'))
Observable('obs_EGFRpb', EGFR(k='p',ec='b'))
Observable('obs_EGFRp', EGFR(k='p',ec='ub'))
Observable('obs_sEGFR', EGFR())
