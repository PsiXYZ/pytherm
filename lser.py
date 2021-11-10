import db.lser.cons as cons
import db.lser.io as io
import sm


def calc_desc():
    pass


def calc_volume(formula: str, bound_count: int) -> float:
    volume = 0
    el, mat = sm.get_el([formula])
    for i in range(len(el)):
        volume += cons.atomic_volumes[el[i]] * mat[i][0]
    return (volume - cons.bond_const*bound_count)/100


def calc_all_subs():
    subs = io.read_subs()
    data = [0] * len(subs)
    for i in range(len(subs)):
        gr = {}
        buf = [0] * 7
        E,S,A,B,V,L = 0,0,0,0,0,0

        # V
        V = calc_volume(subs[i][1], int(subs[i][2]))

        buf = subs[i][3].split(" ")
        buf2 = []
        for j in buf:
            buf2 = j.split('*')
            gr[buf2[1]] = int(buf2[0])

        for j in gr:
            E += cons.c1[j][0] * gr[j]
            S += cons.c1[j][1] * gr[j]
            B += cons.c1[j][2] * gr[j]
            L += cons.c1[j][3] * gr[j]
        E += cons.c1['intercept'][0]
        S += cons.c1['intercept'][1]
        B += cons.c1['intercept'][2]
        L += cons.c1['intercept'][3]

        if int(subs[i][4]) != 0:
            gr = {}
            buf = subs[i][4].split(" ")
            for j in buf:
                buf2 = j.split('*')
                gr[buf2[1]] = buf2[0]
            for j in gr:
                A += cons.c2[j] * gr[j]
            A += cons.c2['intercept']
        data[i] = ([subs[i][0], E, S, A, B, V, L])

    io.write_subs(data)

calc_all_subs()