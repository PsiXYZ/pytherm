from pytherm.activity import unifac as uf
from datetime import datetime
import numpy as np

n = 1_000

params = uf.datasets.get_DOR()
c = (0.5, 0.5)
subs = {
    "hexane": "2*CH3 4*CH2",
    "ethanol": "1*CH3 1*CH2 1*OH(P)",
}
s = uf.SubstancesUNIFAC()
s.get_from_dict(subs)
am = uf.UNIFAC(params, s)
am.get_y(c, 298)
start_time = datetime.now()
for i in range(n):
    am.get_y(c, 298)
print(datetime.now() - start_time)


from pytherm.activity import unifac_numba as uf

substances = uf.datasets.SubstancesUNIFAC()
substances.get_from_dict(subs)
am = uf.UNIFAC(dataset=uf.datasets.DOR, substances=substances)
c = np.array(c)
am.get_y(c, 298)
start_time = datetime.now()
for i in range(n):
    am.get_y(c, 298)
print(datetime.now() - start_time)

from pytherm.activity import unifac_old as uf

substances = uf.datasets.SubstancesUNIFAC()
substances.get_from_dict(subs)
am = uf.UNIFAC(dataset=uf.datasets.DOR, substances=substances, dict_mode=False)
am.get_y(c, 298)
start_time = datetime.now()
for i in range(n):
    am.get_y(c, 298)
print(datetime.now() - start_time)