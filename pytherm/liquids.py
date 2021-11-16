from math import log, log10
from .solver import minimize


def get_yinf(am, ph, comp_name):
    n_ph = 1
    n_comp = 1E-9
    x = {}
    n = {}
    for i in ph:
        n[i] = ph[i] * n_ph
    n[comp_name] = n_comp
    for i in ph:
        x[i] = n[i] / n_ph

    # ph[comp_name] = 1E-10
    y = am.get_y(x)
    # print(x, "--\n",y)

    return y[comp_name]


def get_kexp(am, ph1, ph2, comp_name):
    y1 = get_yinf(am, ph1, comp_name)
    y2 = get_yinf(am, ph2, comp_name)
    return y2 / y1


def get_kp(am, ph1, ph2, r, comp_name):

    v1 = 0
    v2 = 0
    for i in ph1.keys():
        if i != comp_name:
            v1 += ph1[i] / r[i]
            v2 += ph2[i] / r[i]
    kexp = get_kexp(am, ph1, ph2, comp_name)
    rez = kexp * v2/v1
    print("V2/V1 = ", v2/v1)

    return rez


def lle_notifier(f, ph1, ph2):
    comp = list(ph1.keys())
    n_comp = len(comp)
    # ВЫВОД
    print("F", round(log10(f)))
    n_vals = 4
    for s in comp:
        print("{} {:<12.12s} {:<7.4f} {:<7.4f}".format("x", s, round(ph1[s], n_vals), round(ph2[s], n_vals)))

    s1, s2 = 0, 0
    for i in range(n_comp):
        s1 += ph1[comp[i]]
        s2 += ph2[comp[i]]
    print("{:14} {:<7.4f} {:<7.4f}".format("s", round(s1, n_vals), round(s2, n_vals)))


def find_lle(ph1, ph2, am, notifier=lle_notifier):
    # ph1 -> ph2
    n1 = 10
    n2 = 10
    min_ksi = 1e-08  # минимально возможные начения ksi

    comp = list(ph1.keys())
    n_comp = len(comp)

    bound = []
    # расчет n0
    n0_1 = [0] * n_comp
    n0_2 = [0] * n_comp
    for i in range(n_comp):
        if ph1[comp[i]] == 0:
            n0_1[i] = 0
            n0_2[i] = ph2[comp[i]] * n2
            bound.append([- ph2[comp[i]] * n2 + min_ksi, - min_ksi])
        else:
            n0_1[i] = ph1[comp[i]] * n1
            n0_2[i] = ph2[comp[i]] * n2
            bound.append([min_ksi, ph1[comp[i]] * n1 - min_ksi])

    r_mat = []
    for i in range(n_comp):
        r_mat.append([0] * n_comp)
    for i in range(n_comp):
        r_mat[i][i] = 1

    def get_loga(k):
        n1 = 0
        n2 = 0
        n_1 = [0] * n_comp
        n_2 = [0] * n_comp
        a1 = [0] * n_comp
        a2 = [0] * n_comp
        # расчет n при ksi
        for i in range(n_comp):
            n_1[i] = n0_1[i] - k[i]
            n_2[i] = n0_2[i] + k[i]
        for i in range(n_comp):
            n1 += n_1[i]
            n2 += n_2[i]
        # расчет х
        for i in range(n_comp):
            ph1[comp[i]] = n_1[i] / n1
            ph2[comp[i]] = n_2[i] / n2
        y1 = am.get_y(ph1)
        y2 = am.get_y(ph2)
        for i in range(n_comp):
            a1[i] = ph1[comp[i]] * y1[comp[i]]
            a2[i] = ph2[comp[i]] * y2[comp[i]]

        lg = [0] * n_comp
        for i in range(n_comp):
            lg[i] = log(a2[i] / a1[i])
        return lg

    f = minimize(get_loga, bound)[1]

    if notifier is not None:
        lle_notifier(f, ph1, ph2)


def calculate_kp_fromsubs(ph1, ph2, subs):
    pass


'''
def fin(ph1, ph2, am):
    # ph1 -> ph2
    n1 = 10
    n2 = 10
    cer = 1e-08
    fabs = 1e-5

    ksi = [0] * len(ph1)

    comp = list(ph1.keys())
    n_comp = len(comp)

    #двигаю кси для начала оптимизации
    ksi = [0] * n_comp


    # расчет n0
    n0_1 = [0] * n_comp
    n0_2 = [0] * n_comp
    for i in range(n_comp):
        if ph1[comp[i]] == 0:
            n0_1[i] = cer
            n0_2[i] = ph2[comp[i]] * n2 - cer
        else:
            n0_1[i] = ph1[comp[i]] * n1 - cer
            n0_2[i] = ph2[comp[i]] * n2 + cer

    # расчет n при начальных ksi
    n_1 = [0] * n_comp
    n_2 = [0] * n_comp
    # расчет n при ksi
    for i in range(n_comp):
        n_1[i] = n0_1[i]
        n_2[i] = n0_2[i]

    # расчет х
    for i in range(n_comp):
        ph1[comp[i]] = n_1[i] / n1
        ph2[comp[i]] = n_2[i] / n2

    while (1):
        a1 = [0] * n_comp
        a2 = [0] * n_comp
        y1 = am.get_y(ph1)
        y2 = am.get_y(ph2)
        for i in range(n_comp):
            a1[i] = ph1[comp[i]] * y1[comp[i]]
            a2[i] = ph2[comp[i]] * y2[comp[i]]
        lg = [0] * n_comp
        for i in range(n_comp):
            lg[i] = abs(log(a2[i] / a1[i]))
        max_i = 0
        for i in range(n_comp):
            if lg[i] > lg[max_i]:
                max_i = i

        cords = [0 - n0_2[max_i], n0_1[max_i]]
        vals = [0, 0]
        while (1):
            n1 = 0
            n2 = 0

            # расчет n при ksi
            for i in range(n_comp):
                n_1[i] = n0_1[i] - ksi[i]
                n_2[i] = n0_2[i] + ksi[i]

            for i in range(n_comp):
                n1 += n_1[i]
                n2 += n_2[i]

            # расчет х
            for i in range(n_comp):
                ph1[comp[i]] = n_1[i] / n1
                ph2[comp[i]] = n_2[i] / n2

            y1 = am.get_y(ph1)
            y2 = am.get_y(ph2)
            a1 = y1[comp[max_i]] * ph1[comp[max_i]]
            a2 = y2[comp[max_i]] * ph2[comp[max_i]]
            # print(lg, max_i + 1, a1, a2)
            print(ksi, log(a2 / a1))
            #если >0 то нужно уменьшать
            if log(a2 / a1) > 0:
                buf = ksi[max_i]
                ksi[max_i] = (cords[0] + ksi[max_i]) / 2
                cords[1] = buf
            else:
                buf = ksi[max_i]
                ksi[max_i] = (ksi[max_i] + cords[1]) / 2
                cords[0] = buf
            if abs(log(a2 / a1)) < fabs/10:
                break

        if sum(lg) < fabs:
            break

    # ВЫВОД
    print("F", round(log10(sum(lg))))
    n_vals = 4
    for s in comp:
        print("{} {:<12.12s} {:<7.4f} {:<7.4f}".format("x", s, round(ph1[s], n_vals), round(ph2[s], n_vals)))

    s1, s2 = 0, 0
    for i in range(n_comp):
        s1 += ph1[comp[i]]
        s2 += ph2[comp[i]]
    print("{:14} {:<7.4f} {:<7.4f}".format("s", round(s1, n_vals), round(s2, n_vals)))
'''