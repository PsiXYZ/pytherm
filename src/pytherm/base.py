from math import log, exp
from .constants import R


def get_Ht_kelly(h0, kelly, T):
    h = (h0 * 1000
         + kelly[0] * (T - 298)
         + kelly[1] * (T ** 2 - 298 ** 2) / 2 * 10 ** -3
         + kelly[2] * (T ** 3 - 298 ** 3) / 3 * 10 ** -6
         - kelly[3] * (1 / T - 1 / 298) * 10 ** 5)
    return h


def get_St_kelly(s0, kelly, T):
    s = (s0
         + kelly[0] * log(T / 298)
         + kelly[1] * (T - 298) * 10 ** -3 + kelly[2] * (T ** 2 - 298 ** 2) / 2 * 10 ** -6
         - kelly[3] * (1 / (T ** 2) - 1 / (298 ** 2)) / 2 * 10 ** 5)
    return s


def get_Gt_kelly(h0, s0, kelly, T):
    H = get_Ht_kelly(h0, kelly, T)
    S = get_St_kelly(s0, kelly, T)
    return H - T * S


def get_k(G, T, r_mat):
    k = []
    for r_v in r_mat:
        s = 0
        for i in range(len(r_v)):
            s += r_v[i] * G[i]
        k.append(exp(-s / (R * T)))
    return k


def gess(r_v, fun):
    s = 0
    for i in range(len(r_v)):
        s += r_v[i] * fun[i]
    return s
