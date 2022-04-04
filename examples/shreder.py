from math import e
import pytherm.activity.unifac as uf
import numpy as np

Mw = 202
T_m = 353.8
T = 298
R = 8.314
H_m = 19000
lnx = H_m  / R * (1/T_m  - 1/T)

comp = ['hexane', 'naphthalene']
mode = "DOR"
am = uf.Unifac(comp, unifac_mode=mode)
ph1 = {}
ph1[comp[0]] = 1
ph1[comp[1]] = 0

conc = np.linspace(0.01, 0.99, 1000)
mse = 10000
sol = 0
for x in conc: 
    ph1[comp[1]] = x
    ph1[comp[0]] = 1 - x
    y = am.get_y(ph1, temperature=T)
    if (ph1[comp[1]] * y[comp[1]] - e ** lnx) ** 2 < mse:
        mse = (ph1[comp[1]] * y[comp[1]] - e ** lnx) ** 2
        sol = ph1[comp[1]]
    print(ph1[comp[1]], ph1[comp[1]] * y[comp[1]], e ** lnx)
print(sol)