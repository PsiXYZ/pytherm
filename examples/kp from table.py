import pytherm.lle as lq
import pytherm.activity.unifac_old as uf
import csv

inp_file = "unifac in2 t.csv"
out_file = "unifac out t.csv"
comp = ["hexane", "acetonitrile"]
mode_unifac = "INF"
mode_eq = False # True для расчета равновесного состава по юнифаку

ph1 = {}
ph2 = {}
ph1[comp[0]] = 0.945
ph1[comp[1]] = 1 - ph1[comp[0]]
ph2[comp[0]] = 0.050
ph2[comp[1]] = 1 - ph2[comp[0]]
r = {}
r[comp[0]] = 0.6548 / 86
r[comp[1]] = 0.786 / 41


with open(inp_file, 'r') as df:
    reader = csv.reader(df, delimiter = ';', lineterminator= '\n')
    subs = []
    for row in reader:
        subs.append(row[0])

if mode_eq:
    am = uf.Unifac(comp, unifac_mode=mode_unifac)
    ph1 = {}
    ph2 = {}

    ph1['hexane'] = 1
    ph1['acetonitrile'] = 0

    ph2['hexane'] = 0
    ph2['acetonitrile'] = 1

    lq.find_lle(ph1, ph2, am)


out_kp = {}
for s in subs:
    p1 = ph1.copy()
    p2 = ph2.copy()

    p1[s] = 0
    p2[s] = 0


    am = uf.Unifac(p1, unifac_mode=mode_unifac)
    rez = lq.get_kp(am, p1, p2, r, s)
    out_kp[s] = rez

with open(out_file, 'w') as df:
    writer = csv.writer(df, delimiter = ';', lineterminator= '\n')
    writer.writerow([mode_unifac, "EQ CALC" if mode_eq else "EX DATA", ""])
    writer.writerow(['eq', 'x', ""])
    writer.writerow([comp[0],ph1[comp[0]], ph1[comp[1]]])
    writer.writerow([comp[1],ph2[comp[0]], ph2[comp[1]]])

    writer.writerow(['sub', f'{mode_unifac}, {mode_eq}', ""])
    for i in subs:
        writer.writerow([i, out_kp[i], ""])