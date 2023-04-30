import pytherm.activity.unifac as uf

phase = {
    'hexane': 0.99,
    'acetonitrile': 0.01
}
subs_dict = {
    "hexane": "2*CH3 4*CH2",
    "acetonitrile": "1*CH3CN",
}

subs = uf.datasets.SubstancesUNIFAC()
subs.get_from_dict(subs_dict)
am = uf.UNIFAC(dataset=uf.datasets.DOR, substances=subs, dict_mode=True)
print(am.get_y(phase))
