import pytherm.lle as lq
import pytherm.activity.unifac as uf


comp = ["hexane", "acetonitrile"]
target = "naphthalene"

r = {}  # density/Mw
r["hexane"] = 0.6548 / 86
r["acetonitrile"] = 0.786 / 41

subs_dict = {
    "hexane": "2*CH3 4*CH2",
    "acetonitrile": "1*CH3CN",
    "naphthalene": "8*ACH 2*AC"
}

ph1 = {}
ph1["hexane"] = 0.9532
ph1["acetonitrile"] = 1 - ph1["hexane"]
ph1[target] = 0

ph2 = {}
ph2["hexane"] = 0.0629
ph2["acetonitrile"] = 1 - ph2["hexane"]
ph2[target] = 0


subs = uf.datasets.SubstancesUNIFAC()
subs.get_from_dict(subs_dict)
am = uf.UNIFAC(dataset=uf.datasets.VLE, substances=subs, dict_mode=True)
rez = lq.get_kp(am, ph1, ph2, r, target)
print(rez)
