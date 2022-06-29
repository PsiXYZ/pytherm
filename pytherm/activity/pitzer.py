from math import exp, sqrt
import numpy as np
from ..sm import get_charge_dict
from .db import pitzer as datasets


class Pitzer:
    ph = []
    cations = []
    anions = []
    neutral = []
    charge = {}

    db: datasets.ParametersPitzer

    def __init__(self, ph, db):
        self.charge = get_charge_dict()
        self.db = db

        self.ph = list(ph.keys())
        for s in ph:
            if self.charge[s] > 0:
                self.cations.append(s)
            elif self.charge[s] < 0:
                self.anions.append(s)
            else:
                self.neutral.append(s)

    def get_f(self, I, A, b=1.2):
        return - 4 * A * I * np.log(1 + b * np.sqrt(I)) / b

    def get_A(self, ro=0.99707, D=78.4, T=298):
        return (0.13422 *
                (4.1725332 - 0.1481291 * T ** (0.5)
                 + 1.5188505 * 10 ** (-5) * T ** 2
                 - 1.8016317 * 10 ** (-8) * T ** 3
                 + 9.3816144 * 10 ** (-10) * T ** (3.5)))

    def get_I(self, ph):
        I = 0
        for s in ph:
            I += ph[s] * self.charge[s] ** 2
        return I / 2

    def get_B(self, s1, s2, I):
        if self.db['ca'][s1][s2][-1] == 0:
            b0 = self.db['ca'][s1][s2][0]
            b1 = self.db['ca'][s1][s2][1]
            a1 = self.db['ca'][s1][s2][4]
            g1 = self.get_g(a1 * np.sqrt(I))
            return b0 + b1 * g1
        else:
            b0 = self.db['ca'][s1][s2][0]
            b1 = self.db['ca'][s1][s2][1]
            b2 = self.db['ca'][s1][s2][2]
            a1 = self.db['ca'][s1][s2][4]
            a2 = self.db['ca'][s1][s2][5]

            g1 = self.get_g(a1 * np.sqrt(I))
            g2 = self.get_g(a2 * np.sqrt(I))

            return b0 + b1 * g1 + b2 * g2

    def get_g(self, x):
        return 2 * (1 - (1 + x) * exp(-x)) / x ** 2

    def get_C(self, s1, s2):
        C_f = self.db['ca'][s1][s2][3]
        z1 = self.charge[s1]
        z2 = self.charge[s2]

        return C_f / (2 * sqrt(abs(z1 * z2)))

    def get_Z(self, ph):
        z = 0
        for s in ph:
            z += np.abs(self.charge[s]) * ph[s]
        return z

    def get_gibbs(self, ph):
        def get_ca():
            value_ca = 0
            for c in self.cations:
                m_c = ph[c]
                for a in self.anions:
                    if self.db.is_exist('ca', c, a):
                        m_a = ph[a]
                        B = self.get_B(c, a, I)
                        C = self.get_C(c, a)
                        value_ca += m_c * m_a * (2 * B + Z * C)
            return value_ca

        def get_cc():
            value_cc = 0
            for c1 in self.cations:
                m_c1 = ph[c1]
                for c2 in self.cations:
                    if self.db.is_exist('cc', c1, c2):
                        m_c2 = ph[c2]
                        fi = self.db['cc'][c1][c2] + \
                            self.get_theta_e(c1, c2, I)
                        value_cc += m_c1 * m_c2 * fi
            return value_cc

        def get_aa():
            value_aa = 0
            for a1 in self.cations:
                m_a1 = ph[a1]
                for a2 in self.cations:
                    if self.db.is_exist('aa', a1, a2):
                        m_a2 = ph[a2]
                        fi = self.db['aa'][a1][a2] + \
                            self.get_theta_e(a1, a2, I)
                        value_aa += m_a1 * m_a2 * fi
            return value_aa

        def get_cca():
            value_cca = 0
            for c1 in self.cations:
                m_c1 = ph[c1]
                for c2 in self.cations:
                    m_c2 = ph[c2]
                    for a in self.anions:
                        if self.db.is_exist('cca', 'c1', 'c2', 'a'):
                            m_a = self.ph[a]
                            psi = self.db['cca'][c1][c2][a]
                            value_cca += 0.5 * m_c1 * m_c2 * m_a * psi
            return value_cca

        def get_caa():
            value_caa = 0
            for a1 in self.cations:
                m_a1 = ph[a1]
                for a2 in self.cations:
                    m_a2 = ph[a2]
                    for c in self.anions:
                        if self.db.is_exist('caa', c, a1, a2):
                            m_c = ph[c]
                            psi = self.db['caa'][c][a1][a2]
                            value_caa += 0.5 * m_a1 * m_a2 * m_c * psi
            return value_caa

        I = self.get_I(ph)
        A = self.get_A()
        Z = self.get_Z(ph)
        gibbs = self.get_f(I, A)
        gibbs += get_ca()
        if len(self.cations) > 1:
            gibbs += get_cc()
            gibbs += get_cca()
        if len(self.anions) > 1:
            gibbs += get_aa()
            gibbs += get_caa()
        return gibbs

    def get_harvie_j(self, x):
        akI = (
            -0.000000000010991,
            -0.000000000002563,
            0.000000000001943,
            0.000000000046333,
            -0.000000000050847,
            -0.000000000821969,
            0.000000001229405,
            0.000000013522610,
            -0.000000025267769,
            -0.000000202099617,
            0.000000396566462,
            0.000002937706971,
            -0.000004537895710,
            -0.000045036975204,
            0.000036583601823,
            0.000636874599598,
            0.000388260636404,
            -0.007299499690937,
            -0.029779077456514,
            -0.060076477753119,
            1.925154014814667
        )
        # Values from Table B-1, final column (akII)
        akII = (
            0.000000000237816,
            -0.000000002849257,
            -0.000000006944757,
            0.000000004558555,
            0.000000080779570,
            0.000000216991779,
            -0.000000250453880,
            -0.000003548684306,
            -0.000004583768938,
            0.000034682122751,
            0.000087294451594,
            -0.000242107641309,
            -0.000887171310131,
            0.001130378079086,
            0.006519840398744,
            -0.001668087945272,
            -0.036552745910311,
            -0.028796057604906,
            0.150044637187895,
            0.462762985338493,
            0.628023320520852,
        )
        x_vec = np.full_like(akI, x)
        ak = np.where(x_vec < 1, akI, akII)
        z = np.where(
            x < 1, 4 * x ** 0.2 - 2,
            40 / 9 * x ** -0.1 - 22 / 9
        )
        dz_dx = np.where(
            x < 1, 4 * x ** -0.8 / 5,
            -4 * x ** -1.1 / 9
        )
        b2, b1, b0 = 0.0, 0.0, 0.0
        d2, d1, d0 = 0.0, 0.0, 0.0
        for a in ak:
            b2 = b1 * 1
            b1 = b0 * 1
            d2 = d1 * 1
            d1 = d0 * 1
            b0 = z * b1 - b2 + a  # Eq. (B-23/27)
            d0 = b1 + z * d1 - d2  # Eq. (B-24/28)
        J = 0.25 * x - 1 + 0.5 * (b0 - b2)  # Eq. (B-29)
        Jp = 0.25 + 0.5 * dz_dx * (d0 - d2)  # Eq. (B-30)
        return J, Jp

    def get_x_ij(self, z1, z2, I):
        A = self.get_A()
        return 6 * z1 * z2 * A * np.sqrt(I)

    def get_theta_e(self, s1, s2, I):
        z_m = self.charge[s1]
        z_n = self.charge[s2]

        x_mn = self.get_x_ij(z_m, z_n, I)
        x_mm = self.get_x_ij(z_m, z_m, I)
        x_nn = self.get_x_ij(z_n, z_n, I)
        J1, Jp1 = self.get_harvie_j(x_mn)
        J2, Jp2 = self.get_harvie_j(x_mm)
        J3, Jp3 = self.get_harvie_j(x_nn)

        E = (z_m * z_n / (4 * I)) * (J1 - 0.5 * J2 - 0.5 * J3)
        # Ep = (- E / I + (z_m * z_n / (8 * I ** 2))
        #       * (x_mn * Jp1 - 0.5 * x_mm * Jp2 - 0.5 * x_nn * Jp3))
        return E

    def get_y(self, ph):
        y = {}
        for s in self.ph:
            y[s] = np.exp(self.grad(ph, s))
        return y

    def grad(self, ph, s, dm=1e-10):
        ph_b = {}
        for i in self.ph:
            ph_b[i] = ph[i]

        ph_b[s] = ph[s] - dm
        y1 = self.get_gibbs(ph_b)
        ph_b[s] = ph[s] + dm
        y2 = self.get_gibbs(ph_b)

        return (y2 - y1) / (2 * dm)

    def get_a_water(self, ph, n_m=55.50837):
        osmotic = self.get_osmotic()
        s_m = 0
        for i in ph:
            s_m += ph[i]
        lna = - osmotic / (n_m / s_m)
        return np.exp(lna)

    def get_osmotic(self, ph, dw=1e-4, R=8.3145, T=298):
        ph1 = {}
        ph2 = {}
        w_water = 1

        n_i = {}
        for i in ph:
            n_i[i] = ph[i]

        for s in ph:
            ph1[s] = n_i[s] / (w_water - dw)
        self.ph = ph1
        y1 = self.get_gibbs() * (w_water - dw) * R * T

        for s in ph:
            ph2[s] = n_i[s] / (w_water + dw)
        self.ph = ph2
        y2 = self.get_gibbs() * (w_water + dw) * R * T

        self.ph = ph
        dG_dw = (y2 - y1) / (2 * dw)
        s_m = 0
        for i in ph:
            s_m += ph[i]
        return 1 - dG_dw / (R * T * s_m)


# m_nacl = 3.9
# m_naoh = 4
# ph13 = {
#     "Na": m_nacl + m_naoh,
#     "Cl": m_nacl,
#     "OH": m_naoh,
# }

# am = Pitzer(ph13)
# y = am.get_y(ph13)
# print(y)
# # print(
# #     ph['Na'] * ph['Cl'] * y['Na'] * y['Cl']
# # )


# def set_ph(m_nacl, m_naoh=0):
#     ph = {
#         "Na": m_nacl + m_naoh,
#         "Cl": m_nacl,
#         "OH": m_naoh,
#     }
#     return ph


# def PR(ph, y):
#     return ph['Na'] * ph['Cl'] * y['Na'] * y['Cl']


# def f(x):
#     ph = set_ph(x)
#     y = am.get_y(ph)
#     x = PR(ph, y)
#     return x - 37.18


# def opti(bounds=(0.1, 8),
#          f_tol=1e-18):
#     a = bounds[0]
#     b = bounds[1]
#     while(1):
#         x = (a + b) / 2
#         f_a = f(a)
#         f_b = f(b)
#         f_x = f(x)
#         if f_x * f_a < 0:
#             b = x
#         elif f_x * f_b < 0:
#             a = x
#         print(x, a, b, f(x))
#         if abs(b - a) < f_tol or abs(f(x)) < 1e-2:
#             print(x)
#             break


# # opti()
