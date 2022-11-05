from math import log


def find_eq(r_mat, n0, K, T, P, ftol=1e-12, fabs=1e-4, notifier=None):
    # транспонирование r_mat в а
    a = []
    for i in range(len(r_mat[0])):
        b = []
        for j in range(len(r_mat)):
            b.append(r_mat[j][i])
        a.append(b)

    # пересчет к
    for i in range(len(K)):
        K[i] = K[i] * P ** sum(r_mat[i])

    def get_ns(ksi):
        n = []
        for i in range(len(a)):
            s = 0
            for j in range(len(ksi)):
                s += ksi[j] * a[i][j]
            n.append(n0[i] + s)

        return n

    def get_pr(ksi):
        pr = []
        xi = []
        ns = get_ns(ksi)
        s = sum(ns)
        for i in ns:
            xi.append(i / s)
        for i in r_mat:
            s1 = 1
            for j in range(len(i)):
                s1 *= xi[j] ** i[j]
            pr.append(s1 * P ** sum(i))
        return pr

    def fi(ksi):
        res = []
        pr = get_pr(ksi)
        for i in range(len(K)):
            res.append(abs(8.314 * T * (log(pr[i]) - log(K[i]))))
        return res

    def uf(ksi, r):
        n = get_ns(ksi)
        xi = []
        s = sum(n)
        for i in n:
            xi.append(i / s)

        x = 1
        y = 1
        for i in range(len(r_mat[r])):
            if r_mat[r][i] > 0:
                x *= xi[i] ** r_mat[r][i]
            if r_mat[r][i] < 0:
                y *= xi[i] ** (- r_mat[r][i])
        return y * K[r] - x

    def get_rl(ri):
        p_ksi[ri] = 0
        n = get_ns(p_ksi)
        vals = []
        for i in range(len(n0)):
            if r_mat[ri][i] > 0:
                buf = n[i] / r_mat[ri][i]
                vals.append(- buf)
        r_l = max(vals)
        return r_l

    def get_rr(ri):
        p_ksi[ri] = 0
        n = get_ns(p_ksi)
        vals = []
        for i in range(len(n0)):
            if r_mat[ri][i] < 0:
                buf = n[i] / r_mat[ri][i]
                vals.append(- buf)
        r_r = min(vals)
        return r_r

    ksi = [0] * len(K)
    while (1):
        f = fi(ksi)
        print("F = ", f)
        print("ksi = ", ksi)

        if sum(f) < fabs:
            break

        v = -1
        ri = 0
        # поиск реакции с мах f
        for i in range(len(f)):
            if f[i] > v:
                v = f[i]
                ri = i

        # расчет n при ksi[i] = 0
        p_ksi = ksi
        p_ksi[ri] = 0
        n = get_ns(p_ksi)

        r_l = get_rl(ri)
        r_r = get_rr(ri)

        rs = [r_l, 0, r_r]

        # начало оптимизации
        vals_u = [0, 0]
        while (1):
            rs[1] = (rs[0] + rs[2]) / 2
            p_ksi[ri] = rs[1]
            vals_u[1] = vals_u[0]
            vals_u[0] = uf(p_ksi, ri)

            if vals_u[0] < 0:
                rs[2] = rs[1]
            else:
                rs[0] = rs[1]

            if vals_u[0] == 0:
                ksi[ri] = rs[1]
                break
            if abs((vals_u[0] - vals_u[1]) / vals_u[0]) < ftol:
                ksi[ri] = rs[1]
                break

    print("final = ", ksi)
    print("Fs final = ", fi(ksi))
    print("abs F final = ", sum(fi(ksi)))
    print("PR final = ", get_pr(ksi))
    print("Ns", get_ns(ksi))
