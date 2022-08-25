import pytherm.activity.unifac_old as uf
from math import log

T_m = 360.4
T = 298
R = 8.314
H_m = 23550
lnx = H_m / R * (1/T_m - 1/T)

comp = ['water', 'DDT']
mode = "VLE"
am = uf.Unifac(comp, unifac_mode=mode)
ph1 = {}
ph1[comp[0]] = 1
ph1[comp[1]] = 0

bounds = (1e-20, 0.99)
f_tol = 1e-18
a = bounds[0]
b = bounds[1]


def f(x):
    ph1[comp[1]] = x
    ph1[comp[0]] = 1 - x
    y = am.get_y(ph1, temperature=T)
    return (log(ph1[comp[1]] * y[comp[1]]) - lnx)


while(1):
    x = (a + b) / 2
    if f(x) * f(a) < 0:
        b = x
    elif f(x) * f(b) < 0:
        a = x
    # print(x, a, b, f(x))
    if abs(b - a) < f_tol or abs(f(x)) < 1e-10:
        break

print(f"X: {x:.1e} \tM: {x * 354 / 18:.1e} \tD: {f(x):.1e}")
# print("x", ph1)
# print("y", am.get_y(ph1))
