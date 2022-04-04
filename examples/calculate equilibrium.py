import pytherm.liquids as lq
import pytherm.activity.unifac as uf

comp = ['hexane', 'acetonitrile']
target = "toluene"
mode = "VLE"
am = uf.Unifac(comp, unifac_mode=mode)

ph1 = {}
ph2 = {}

ph1['hexane'] = 1
ph1['acetonitrile'] = 0

ph2['hexane'] = 0
ph2['acetonitrile'] = 1

lq.find_lle(ph1, ph2, am, temperature=298)
