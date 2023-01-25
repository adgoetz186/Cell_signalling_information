from pysb import *

# Parameters in order:
# kprod - unbound unphosphorylated receptor synthesis rate
# kdeg - all receptors degredation rate
# kbind - receptor ligand binding rate constant
# kunbind - receptor ligand unbinding rate constant
# krp - receptor phosphorylation rate constant
# krdp - receptor dephosphorylation rate constant
# kadp - Akt dephosphorylation rate constant
# kap - Akt phosphorylation rate constant
# kfp - FoxO phosphorylation rate constant
# kfdp - FoxO dephosphorylation rate constant
# kinflux - FoxO influx into nucleus rate constant
# kefflux - FoxO efflux out of nucleus rate constant
# IGFR_0 - Initial unphosphorylated unbound IGFR count
# Akt_0 - Initial unphosphorylated Akt count
# FoxO_c_0 - Initial unphosphorylated FoxO in the cytoplasm
# FoxO_n_0 - Initial unphosphorylated FoxO in the nucleus
# IGF - Ligand concentration UNITS MATTER, value MUST be in nM

Model()
Compartment('e', dimension=2, parent=None)
Compartment('c',dimension=3,parent=e)
Compartment('ne',dimension=2,parent=c)
Compartment('n',dimension=3,parent=ne)

# Due to the replacement rate of the ligand, the ligand concentration acts like a constant
# so IGF binding acts as a change of state of IGFR rather than a binding event
Monomer('IGFR', ['k','ec'],{'k':['up','p'],'ec':['b','ub']})
Monomer('Akt', ['k'],{'k':['up','p']})
Monomer('FoxO', ['k'],{'k':['up','p']})

Parameter('kprod',1.0)
Parameter('kdeg',1.0)

Parameter('kbind',1.0)
Parameter('kunbind',1.0)
Parameter('krp',1.0)
Parameter('krdp',1.0)
Parameter('kadp',1.0)
Parameter('kap',1.0)
Parameter('kfp',1.0)
Parameter('kfdp',1.0)
Parameter('kinflux',1.0)
Parameter('kefflux',1.0)
Parameter('IGFR_0', 1)
Parameter('Akt_0', 1)
Parameter('FoxO_c_0', 1)
Parameter('FoxO_n_0', 1)

Parameter('IGF', 1.0)

# this is effective binding constant, should be igf concentration and rate constant
Expression("keffbind",kbind*IGF)

Initial(IGFR(k='up',ec='ub')**e, IGFR_0)
Initial(Akt(k='up')**c, Akt_0)
Initial(FoxO(k='up')**c, FoxO_c_0)
Initial(FoxO(k='up')**n, FoxO_n_0)

Rule('r_prod', None >> IGFR(k='up',ec='ub')**e,kprod)
Rule('r_deg', IGFR()**e >> None,kdeg)
Rule('r_bind', IGFR(k='up',ec='ub')**e | IGFR(k='up',ec='b')**e,keffbind,kunbind)
Rule('rb_phos', IGFR(k='up',ec='b')**e | IGFR(k='p',ec='b')**e,krp,krdp)
Rule('akt_phos', Akt(k='up')**c + IGFR(k='p',ec='b')**e >> Akt(k='p')**c + IGFR(k='p',ec='b')**e,kap)
Rule('akt_dephos', Akt(k='p')**c >> Akt(k='up')**c,kadp)
Rule('foxo_phos', FoxO(k='up')**c + Akt(k='p') >> FoxO(k='p')**c + Akt(k='p'),kfp)
Rule('foxo_dephos', FoxO(k='p')**c >> FoxO(k='up')**c, kfdp)
Rule('foxo_influx', FoxO(k='up')**c | FoxO(k='up')**n, kinflux,kefflux)

Observable('obs_FoxO_n', FoxO(k='up')**n)
Observable('obs_IGFR', IGFR()**e)
Observable('obs_FoxO_c', FoxO()**c)
