import pytherm.activity.unifac as uf
from math import e, cos, sin, pi
import numpy as np

comp = ["fr1", "fr2", "as2"]
d = {
    "naphthalene": [353.5, 19550, 128.17052],
    "anthracene": [491, 31500, 178.234],
    "phenanthrene": [367.6, 16600, 178.234],
    "pyrene": [422.4, 16700, 202.25],
    "as1": [2510, 17000, 1518],
    "as1/2": [1307, 17000, 759],
    "as1*2": [4903, 17000, 3036],
    "as1-pyr": [2510, 17000, 1518],
    "as1-ovy": [2510, 17000, 1518],
    "as2": [4058, 17000, 2500],
    "as2-pyr": [4100, 17000, 2500],
    "as2-ovy": [4100, 17000, 2500],
    "as3": [827, 17000, 450],
}
target = comp[2]
mode = "VLE"
T = 373  # 298 323 348 373 

Mw = d[target][2]
T_m = d[target][0]
H_m = d[target][1]
R = 8.314
lnx = H_m / R * (1/T_m - 1/T)

am = uf.Unifac([comp[0], comp[2]], unifac_mode=mode)
ph1 = {}
ph1[comp[0]] = 0
ph1[comp[2]] = 0
conc = np.linspace(0.0001, 0.9999, 1000)
mse = 10000
sol = 0
yy = 0
for x in conc: 
    ph1[comp[2]] = x
    ph1[comp[0]] = 1 - x
    y = am.get_y(ph1, temperature=T)
    if (ph1[comp[2]] * y[comp[2]] - e ** lnx) ** 2 < mse:
        mse = (ph1[comp[2]] * y[comp[2]] - e ** lnx) ** 2
        sol = ph1[comp[2]]
        yy = y[comp[2]]
print(sol)

am = uf.Unifac(comp, unifac_mode=mode)
ph = {}
ph[comp[0]] = 1 - sol
ph[comp[1]] = 0
ph[comp[2]] = sol

of = open("sol t.txt", 'w')
of.write(f"{comp[0]}, {comp[1]}, {comp[2]},{comp[0]}, {comp[1]}, {comp[2]}\n")

n = 50
r = 0.01
a = pi/n
v1 = [r, 0] # 2 3
v2 = [v1[0] * cos(-pi/2) - v1[1] * sin(-pi/2), 
        v1[0] * sin(-pi/2) + v1[1] * cos(-pi/2)]
while(1):
    mse = 10000
    x2_min = 0
    x3_min = 0
    ph1 = {}
    ph1[comp[0]] = 0
    ph1[comp[1]] = 0
    ph1[comp[2]] = 0
    for i in range(1, n):
        angle = a * i
        x2_n = ph[comp[1]] + v2[0] * cos(angle) - v2[1] * sin(angle)
        x3_n = ph[comp[2]] + v2[0] * sin(angle) + v2[1] * cos(angle)
        x1_n = 1 - x2_n - x3_n
        
        ph1[comp[0]] = x1_n
        ph1[comp[1]] = x2_n
        ph1[comp[2]] = x3_n

        y = am.get_y(ph1, temperature=T)
        if (ph1[comp[2]] * y[comp[2]] - e ** lnx) ** 2 < mse:
            mse = (ph1[comp[2]] * y[comp[2]] - e ** lnx) ** 2
            x2_min = x2_n
            x3_min = x3_n
            x1_min = 1 - x2_n - x3_n
            yy = y[comp[2]]
    if min(x1_min, x2_min, x3_min) < 0:
        break
    v1 = [x2_min - ph[comp[1]], x3_min - ph[comp[2]]]
    v2 = [v1[0] * cos(-pi/2) - v1[1] * sin(-pi/2), v1[0] * sin(-pi/2) + v1[1] * cos(-pi/2)]
    ph[comp[0]] = x1_min
    ph[comp[1]] = x2_min
    ph[comp[2]] = x3_min
    print(x1_min, x2_min, x3_min, yy * x3_min, e ** lnx, mse)

    sM = x1_min*226 + x2_min*162 + x3_min*Mw
    of.write(f"{x1_min}, {x2_min}, {x3_min},{x1_min*228.6/sM}, {x2_min*161.7/sM}, {x3_min*Mw/sM}\n")
print("________________________DONE________________________")