import pytherm.liquids as lq
import pytherm.activity.unifac as uf

comp = ["hexane", "acetonitrile"]
target = "naphthalene"

r = {}
r["hexane"] = 0.6548 / 86
r["acetonitrile"] = 0.786 / 41

ph1 = {}
ph1["hexane"] = 0.9532
ph1["acetonitrile"] = 1 - ph1["hexane"]
ph1[target] = 0

ph2 = {}
ph2["hexane"] = 0.0629
ph2["acetonitrile"] = 1 - ph2["hexane"]
ph2[target] = 0

am = uf.Unifac(ph1, mode="VLE")
rez = lq.get_kp(am, ph1, ph2, r, target)
print(rez)