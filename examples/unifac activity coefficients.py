import pytherm.activity.unifac as uf

phase = {
    'hexane': 0.99,
    'acetonitrile': 0.01
}
subs_dict = {
    "hexane": "2*CH3 4*CH2",
    "acetonitrile": "1*CH3CN",
}

subs = uf.SubstancesUNIFAC()
subs.get_from_dict(subs_dict)
am = uf.UNIFAC(uf.datasets.DOR(), subs)

print(
    am.get_y(list(phase.values()), 298)
    )
