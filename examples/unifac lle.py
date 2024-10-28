import pytherm.lle as lq
import pytherm.activity.unifac as uf

phase1 = {
    'hexane': 1,
    'acetonitrile': 0
}
phase2 = {
    'hexane': 0,
    'acetonitrile': 1
}
subs_dict = {
    "hexane": "2*CH3 4*CH2",
    "acetonitrile": "1*CH3CN",
}

subs = uf.SubstancesUNIFAC()
subs.get_from_dict(subs_dict)
am = uf.UNIFAC(uf.datasets.DOR(), subs)

lq.find_lle(phase1, phase2, am, T=298)
