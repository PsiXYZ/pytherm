from pytherm.activity import uniquac as uq
T = 25.0 + 273.15
xs = [0.7273, 0.0909, 0.1818]
rs = [0.92, 2.1055, 3.1878]
qs = [1.4, 1.972, 2.4]
inter = [
    [[0, 0], [0, 526.02], [0, 309.64]],
    [[0, -318.06], [0, 0], [0, -91.532]],
    [[0, 1325.1], [0, 302.57], [0, 0]],
]
am = uq.UNIQUAC(rs, qs, inter)
y = am.get_y(xs, T)
print(y)
