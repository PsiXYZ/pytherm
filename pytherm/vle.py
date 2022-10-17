import pytherm.eos as eos
import numpy as np


def fit_Pxy(model: eos.EOS, T: float, points=50):
    system = model.get_system()
    subs = list(system.keys())
    pxy = []
    bubble_curve = [[], []]
    dew_curve = [[], []]

    xi = np.linspace(0.001, 0.999, num=points)
    for x in xi:
        system[subs[0]] = x
        system[subs[1]] = 1 - x
        unst = is_unstable(system, model, T)
        if unst > 0:
            P, rez = find_xy(model=model, T=T, phase_l=system, P_s=unst)
            bubble_curve[0].append(system[subs[0]])
            bubble_curve[1].append(P)

            dew_curve[0].append(rez[subs[0]])
            dew_curve[1].append(P)

            # bubble_curve.append([system[subs[0]], P])
            # dew_curve.append([rez[subs[0]], P])

    return bubble_curve, dew_curve


# def is_unstable(system, model: eos.EOS, T: float, P=1E5):
#     r = model.get_roots(system=system, T=T, P=P)
#     if len(r) == 1:
#         return False
#     else:
#         return True


def is_unstable(system, model: eos.EOS, T: float, P=1E5, points=1000):
    p = []
    b = model.get_b(system)
    vi = np.linspace(b * 1.4, 1600 / 1000000, num=points)
    for i in vi:
        p.append(model.get_p(system=system, T=T, V=i))
    p = np.array(p)
    grad_p = np.gradient(p)
    # i_max = grad_p.argmax()
    if grad_p.max() > 0:
        zeros_index = find_zeros_index(grad_p)
        max_v = p[zeros_index[1]]
        min_v = p[zeros_index[0]]
        if min_v > 0:
            return (max_v + min_v)/2
        else:
            return max_v/2
    else:
        return -1


def find_zeros_index(y):
    zeros = []
    for i in range(1, len(y)):
        if y[i] > 0 and y[i - 1] < 0:
            zeros.append(i)
        if y[i] < 0 and y[i - 1] > 0:
            zeros.append(i)
    return zeros

def find_xy(model: eos.EOS, T: float, phase_l, abs_err=1e-4, P_s=1E5):
    phase_v = {}
    for i in phase_l:
        phase_v[i] = phase_l[i]

    P = P_s
    step = P_s / 4
    err = 2
    prev_err = 3
    while True:
        r = model.get_roots(system=phase_l, T=T, P=P)
        f_l = model.get_f(phase_l, P=P, V=r[0], T=T)

        r = model.get_roots(system=phase_v, T=T, P=P)
        f_v = model.get_f(phase_v, P=P, V=r[-1], T=T)

        yi = {}
        for i in phase_l:
            yi[i] = phase_l[i] * f_l[i] / f_v[i]

        s = 0
        for i in yi:
            s += yi[i]

        for i in yi:
            phase_v[i] = yi[i] / s

        prev_err = err
        err = s - 1
        if abs(err) < abs_err:
            break
        if err / prev_err < 0:
            step /= 2

        if err > 0:
            P += step
        if err < 0:
            if P - step < 0:
                step /= 2
            P -= step

    return P, phase_v
