from datetime import datetime


from pytherm.activity import unifac

params = unifac.datasets.get_DOR()
Mw = (100, 50)
c = (0.5, 0.5)
subs = {
    "hexane": "2*CH3 4*CH2",
    "ethanol": "1*CH3 1*CH2 1*OH(P)",
}
s = unifac.SubstancesUNIFAC()
s.get_from_dict(subs)
am = unifac.UNIFAC_W(params, s, Mw)
y = am.get_y(c, 298)
print(y)


from pytherm.activity import unifac_old

substances = unifac_old.datasets.SubstancesUNIFAC()
substances.get_from_dict(subs)
c = {
    "hexane": 0.5,
    "ethanol": 0.5,
}
Mw = {
    "hexane": 100,
    "ethanol": 50,
}
am = unifac_old.UNIFAC_W(dataset=unifac_old.datasets.DOR, substances=substances, Mw=Mw, dict_mode=True)
y = am.get_y(c, 298)
print(y)