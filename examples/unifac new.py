import pytherm.activity.unifac as uf
uf.datasets.VLE

ph1 = {}
ph1['hexane'] = 0.99
ph1['acetonitrile'] = 1 - ph1['hexane']

subs = uf.datasets.SubstancesUNIFAC()
subs.get_from_defsubs(ph1)

am = uf.Unifac(dataset=uf.datasets.DOR,
               substances=subs)

print(am.get_y(ph1))
