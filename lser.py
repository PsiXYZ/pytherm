import db.lser.constants as cons
import db.lser.io as io
import sm


def calc_desc():
    pass


def calc_volume(formula: str, bound_count: int) -> float:
    volume = 0
    el, mat = sm.get_el([formula])
    for i in range(len(el)):
        volume += cons.atomic_volumes[el[i]] * mat[i][0]
    return (volume - cons.bond_const * bound_count) / 100


def calc_all_subs():
    subs = io.read_subs()
    data = [0] * len(subs)
    for i in range(len(subs)):
        gr = {}
        buf = [0] * 7
        E, S, A, B, V, L = 0, 0, 0, 0, 0, 0

        # V
        V = calc_volume(subs[i][1], int(subs[i][2]))

        buf = subs[i][3].split(" ")
        if buf[-1] == " ":
            buf.remove(buf[-1])
        buf2 = []
        for j in buf:
            buf2 = j.split('*')
            gr[buf2[1]] = int(buf2[0])

        for j in gr:
            E += cons.group_constants[j][0] * gr[j]
            S += cons.group_constants[j][1] * gr[j]
            B += cons.group_constants[j][2] * gr[j]
            L += cons.group_constants[j][4] * gr[j]
        E += cons.group_constants['intercept'][0]
        S += cons.group_constants['intercept'][1]
        B += cons.group_constants['intercept'][2]
        L += cons.group_constants['intercept'][4]

        if not subs[i][4].isalnum():
            gr = {}
            buf = subs[i][4].split(" ")
            for j in buf:
                buf2 = j.split('*')
                gr[buf2[1]] = int(buf2[0])
            for j in gr:
                A += cons.acidity_constants[j] * gr[j]
            A += cons.acidity_constants['intercept']
        data[i] = ([subs[i][0], E, S, A, B, V, L])

    io.write_subs(data)


calc_all_subs()
