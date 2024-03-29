import numpy as np


def find_eq(n0, r_mat, K, T, P, ftol=1e-12, fabs=1e-4, k_lim=10):
    def get_n(ksi):
        n = np.full(len(n0), 0.0)
        for i in range(len(n0)):
            n[i] = r_mat[:, i] @ ksi
        return n0 + n

    def get_pr(ksi):
        pr = np.full(len(K), 0.0)  # реакционные произведения
        ni = get_n(ksi)  # кол-ва вещества i
        n = np.sum(ni)  # общее кол-ва вещества
        xi = ni / n  # мольные доли i

        for i in range(len(K)):
            s = np.power(xi, r_mat[i])
            pr[i] = np.prod(s) * P ** sum(r_mat[i])
        return pr

    def fi(ksi):
        pr = get_pr(ksi)
        # res = (pr - K) ** 2
        res = np.abs(np.log(pr) - np.log(K))
        return res

    def get_rl(p_ksi, ri):
        p_ksi[ri] = 0
        n = get_n(p_ksi)
        vals = []
        for i in range(len(n0)):
            if r_mat[ri][i] > 0:
                buf = n[i] / r_mat[ri][i]
                vals.append(- buf)
        r_l = max(vals)
        return r_l

    def get_rr(p_ksi, ri):
        p_ksi[ri] = 0
        n = get_n(p_ksi)
        vals = []
        for i in range(len(n0)):
            if r_mat[ri][i] < 0:
                buf = n[i] / r_mat[ri][i]
                vals.append(- buf)
        r_r = min(vals)
        return r_r

    # пересчет К на давление
    K = np.array(K)
    for i in range(len(K)):
        K[i] = K[i] * P ** np.sum(r_mat[i])

    # заранее неравновесные реакции проводятся заранее
    for i in range(len(K)):
        if abs(K[i]) < 10 ** (- k_lim):
            p_ksi = np.full(len(K), 0, dtype=float)
            b = get_rl(p_ksi, i)
            p_ksi[i] = b
            n = get_n(p_ksi)
            n0 = n
        if abs(K[i]) > 10 ** k_lim:
            p_ksi = np.full(len(K), 0, dtype=float)
            b = get_rr(p_ksi, i)
            p_ksi[i] = b
            n = get_n(p_ksi)
            n0 = n

    # неравновесные реакции не нужно уравнивать
    # is_optimize = np.full(len(K), True)
    # for i in range(len(K)):
    #     if abs(K[i]) > 10 ** k_lim or abs(K[i]) < 10 ** (- k_lim):
    #         is_optimize[i] = False

    ksi = np.full(len(K), 0, dtype=float)
    while (1):
        f = fi(ksi)
        print("F = ", f)
        print("ksi = ", ksi)

        if np.sum(f) < fabs:
            break

        # поиск реакции с мах f
        v = 0
        ri = 0
        for i in range(len(f)):
            if f[i] > v:
                v = f[i]
                ri = i

        # расчет n при ksi[i] = 0
        p_ksi = ksi
        p_ksi[ri] = 0
        r_l = get_rl(p_ksi, ri)
        r_r = get_rr(p_ksi, ri)
        rs = np.array((r_l, 0, r_r))

        # начало оптимизации
        vals_u = np.array((0.0, 0.0))
        while (1):
            rs[1] = (rs[0] + rs[2]) / 2
            p_ksi[ri] = rs[1]
            # print("F = ", fi(p_ksi))
            # print("Ksi = ", p_ksi)
            res = (K - get_pr(p_ksi))
            vals_u[1] = vals_u[0]
            vals_u[0] = fi(p_ksi)[ri]

            if res[ri] < 0:
                rs[2] = rs[1]
            else:
                rs[0] = rs[1]

            if vals_u[0] == 0:
                ksi[ri] = rs[1]
                break
            if vals_u[0] < fabs / 10:
                ksi[ri] = rs[1]
                break
            tol = np.abs((vals_u[0] - vals_u[1]) / vals_u[0])
            if tol < ftol:
                ksi[ri] = rs[1]
                break

    print("________________DONE__________________________")
    print("KSI =\n", ksi)
    print("F final =\n", fi(ksi))
    print("K =\n", K)
    print("PR final =\n", get_pr(ksi))
    print("Ns =\n", get_n(ksi))
    return get_n(ksi)
