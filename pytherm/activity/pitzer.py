from math import exp, log, sqrt
import db.pitzer.parameters as db
import numpy as np


class Pitzer:
    charge = {}
    cations = []
    anions = []
    neutral = []
    I = 0
    Z = 0
    ph = {}

    pure_substance = {}  # [b0, b1, b2, c, a1, a2]
    theta = {}
    psi = {}
    lamda = {}

    def __init__(self, ph):
        self.charge = db.charge

        for i in db.pure:
            self.pure_substance[f'{i[0]} - {i[1]}'] = i[2:]

        self.theta = db.theta_pp + db.theta_mm
        self.psi = db.psi_ppm + db.psi_mmp
        self.lamda = db.lamda

        self.ph = ph
        for s in ph:
            if self.charge[s] > 0:
                self.cations.append(s)
            elif self.charge[s] < 0:
                self.anions.append(s)
            else:
                self.neutral.append(s)

    def get_A(self, ro=0.99707, D=78.4, T=298):
        # return 1400684*(ro/(D * T)) ** (3/2)
        return 0.13422 * (4.1725332 - 0.1481291 * T ** (0.5) + 1.5188505 * 10 ** (-5) * T ** 2 - 1.8016317 * 10 ** (-8) * T ** 3 + 9.3816144 * 10 ** (-10) * T ** (3.5))

    def get_I(self):
        I = 0
        for s in self.ph:
            I += self.ph[s] * self.charge[s] ** 2
        return I / 2

    def get_F(self):
        def get_f_y():
            I = self.I
            b = 1.2

            f = (sqrt(I) / (1 + b * sqrt(I))) + \
                (2 / b) * log(1 + b * sqrt(I))
            return f

        def get_B_s(s1, s2):
            if self.pure_substance[f'{s1} - {s2}'][-1] == 0:
                b1 = self.pure_substance[f'{s1} - {s2}'][1]
                a1 = self.pure_substance[f'{s1} - {s2}'][4]
                I = self.I
                g1 = get_g_s(a1 * sqrt(I))
                return b1 * g1 / I
            else:
                I = self.I
                b1 = self.pure_substance[f'{s1} - {s2}'][1]
                b2 = self.pure_substance[f'{s1} - {s2}'][2]

                a1 = self.pure_substance[f'{s1} - {s2}'][4]
                a2 = self.pure_substance[f'{s1} - {s2}'][5]

                g1 = get_g_s(a1 * sqrt(I))
                g2 = get_g_s(a2 * sqrt(I))

                return b1 * g1 / I + b2 * g2 / I

        def get_g_s(x):
            return - 2 * (1 - (1 + x + 1/2 * x ** 2) * exp(-x)) / x ** 2

        f_y = get_f_y()
        A = self.get_A()

        ca = 0
        for s1 in self.cations:
            for s2 in self.anions:
                B_s = get_B_s(s1, s2)
                ca += self.ph[s1] * self.ph[s2] * B_s

        cc = 0
        for s1 in self.cations:
            for s2 in self.cations:
                if s1 == s2:
                    continue
                if f'{s1} - {s2}' in self.theta:
                    cc += self.ph[s1] * self.ph[s2] * \
                        self.theta[f'{s1} - {s2}']

        aa = 0
        for s1 in self.anions:
            for s2 in self.anions:
                if s1 == s2:
                    continue
                if f'{s1} - {s2}' in self.theta:
                    aa += self.ph[s1] * self.ph[s2] * \
                        self.theta[f'{s1} - {s2}']

        return - A * f_y + ca + cc + aa

    def get_B(self, s1, s2):
        if self.pure_substance[f'{s1} - {s2}'][-1] == 0:
            I = self.I
            b0 = self.pure_substance[f'{s1} - {s2}'][0]
            b1 = self.pure_substance[f'{s1} - {s2}'][1]
            a1 = self.pure_substance[f'{s1} - {s2}'][4]
            g1 = self.get_g(a1 * sqrt(I))
            return b0 + b1 * g1
        else:
            I = self.I
            b0 = self.pure_substance[f'{s1} - {s2}'][0]
            b1 = self.pure_substance[f'{s1} - {s2}'][1]
            b2 = self.pure_substance[f'{s1} - {s2}'][2]
            a1 = self.pure_substance[f'{s1} - {s2}'][4]
            a2 = self.pure_substance[f'{s1} - {s2}'][5]

            g1 = self.get_g(a1 * sqrt(I))
            g2 = self.get_g(a2 * sqrt(I))

            return b0 + b1 * g1 + b2 * g2

    def get_g(self, x):
        return 2 * (1 - (1 + x) * exp(-x)) / x ** 2

    def get_C(self, s1, s2):
        C_f = self.pure_substance[f'{s1} - {s2}'][3]
        z1 = self.charge[s1]
        z2 = self.charge[s2]

        return C_f / (2 * sqrt(abs(z1 * z2)))

    def get_Z(self):
        z = 0
        for s in self.ph:
            z += abs(self.charge[s]) * self.ph[s]
        return z

    def get_y(self, ph):
        pass

    def get_fi(self, s1, s2):
        if f'{s1} - {s2}' in self.theta:
            return self.theta[f'{s1} - {s2}']
        else:
            return 0

    def get_psi(self, s1, s2, s3):
        if f'{s1} - {s2} - {s3}' in self.psi:
            return self.psi[f'{s1} - {s2} - {s3}']
        else:
            return 0

    def get_lamda(self, s1, s2):
        if f'{s1} - {s2}' in self.lamda:
            return self.lamda[f'{s1} - {s2}']
        else:
            return 0

    def get_lny_m(self, substance):
        self.I = self.get_I()
        self.Z = self.get_Z()
        z_m = self.charge[substance]
        F = self.get_F()

        p1 = 0
        for a in self.anions:
            B = self.get_B(substance, a)
            C = self.get_C(substance, a)
            p1 += self.ph[a] * (2 * B + self.Z * C)

        p2 = 0
        for c in self.cations:
            m = self.ph[c]
            fi = self.get_fi(substance, c)

            buf = 0
            for a in self.anions:
                psi = self.get_psi(substance, c, a)
                buf += self.ph[a] * psi

            p2 += m * (2 * fi + buf)

        p3 = 0
        for a1 in self.anions:
            for a2 in self.anions:
                if a1 == a2:
                    continue
                p3 += self.ph[a1] * self.ph[a2] * \
                    self.get_psi(a1, a2, substance)

        p4 = 0
        for c in self.cations:
            for a in self.anions:
                p4 += self.ph[c] * self.ph[a] * self.get_C(c, a)
        p4 *= abs(z_m)

        return z_m ** 2 * F + p1 + p2 + p3 + p4

    def get_lny_x(self, substance):
        self.I = self.get_I()
        self.Z = self.get_Z()
        z_x = self.charge[substance]
        F = self.get_F()

        p1 = 0
        for c in self.cations:
            B = self.get_B(c, substance)
            C = self.get_C(c, substance)
            p1 += self.ph[c] * (2 * B + self.Z * C)

        p2 = 0
        for a in self.anions:
            m = self.ph[a]
            fi = self.get_fi(substance, a)

            buf = 0
            for c in self.cations:
                psi = self.get_psi(substance, a, c)
                buf += self.ph[a] * psi

            p2 += m * (2 * fi + buf)

        p3 = 0
        for c1 in self.anions:
            for c2 in self.anions:
                if c1 == c2:
                    continue
                p3 += self.ph[c1] * self.ph[c2] * \
                    self.get_psi(c1, c2, substance)

        p4 = 0
        for c in self.cations:
            for a in self.anions:
                p4 += self.ph[c] * self.ph[a] * self.get_C(c, a)
        p4 *= abs(z_x)

        return z_x ** 2 * F + p1 + p2 + p3 + p4

    def get_lny_n(self, substance):
        p1 = 0
        for c in self.cations:
            p1 += self.ph[c] * (2 * self.get_lamda(substance, c))

        p2 = 0
        for a in self.anions:
            p2 += self.ph[a] * (2 * self.get_lamda(substance, a))

    def get_o(self):
        def get_B_f(self, s1, s2):
            I = self.I
            b0 = self.pure_substance[f'{s1} - {s2}'][0]
            b1 = self.pure_substance[f'{s1} - {s2}'][1]
            b2 = self.pure_substance[f'{s1} - {s2}'][2]
            a1 = self.pure_substance[f'{s1} - {s2}'][4]
            a2 = self.pure_substance[f'{s1} - {s2}'][5]

            return b0 + b1 * exp(-a1 * sqrt(I)) + b2 * exp(-a2 * sqrt(I))

        self.I = self.get_I()
        self.Z = self.get_Z()

        p1 = 0
        for i in self.ph:
            p1 += self.ph[i]
        p1 = 2 / p1

        p2 = - (self.get_A() * self.I ** (3 / 2)) / \
            (1 + 1.2 * self.I ** (1 / 2))

        p3 = 0
        for c in self.cations:
            m_c = self.ph[c]
            for a in self.anions:
                B = get_B_f(c, a)
                C = self.get_C(c, a)
                m_a = self.ph[a]
                p3 += m_c * m_a * (B + self.Z * C)

        p4 = 0
        for c1 in self.cations:
            m_c1 = self.ph[c1]
            for c2 in self.cations:
                if c1 == c2:
                    continue
                m_c2 = self.ph[c2]
                fi = self.get_fi(c1, c2)

                buf = 0
                for a in self.anions:
                    m_a = self.ph[a]
                    psi = self.get_psi(c1, c2, a)
                    buf += m_a * psi

                p4 += m_c1 * m_c2 * (fi + buf)

        p5 = 0
        for a1 in self.cations:
            m_a1 = self.ph[a1]
            for a2 in self.cations:
                if a1 == a2:
                    continue
                m_a2 = self.ph[a2]
                fi = self.get_fi(a1, a2)

                buf = 0
                for c in self.anions:
                    m_c = self.ph[c]
                    psi = self.get_psi(a1, a2, c)
                    buf += m_c * psi

                p5 += m_a1 * m_a2 * (fi + buf)

        return p1 + p2 + p3 + p4 + p5

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

    def get_x_ij(self, z1, z2):
        A = self.get_A()
        return 6 * z1 * z2 * A * sqrt(self.I)

    def get_theta_e(self, z_m, z_n):
        x_mn = self.get_x_ij(z_m, z_n)
        x_mm = self.get_x_ij(z_m, z_m)
        x_nn = self.get_x_ij(z_n, z_n)
        J1, Jp1 = self.get_harvie_j(x_mn)
        J2, Jp2 = self.get_harvie_j(x_mm)
        J3, Jp3 = self.get_harvie_j(x_nn)
        I = self.I

        E = (z_m * z_n / (4 * I)) * (J1 - 0.5 * J2 - 0.5 * J3)
        Ep = - E / I + (z_m * z_n / (8 * I ** 2)) * (x_mn * Jp1 - 0.5 * x_mm * Jp2 - 0.5 * x_nn * Jp3)
        return E, Ep


ph = {
    "Mg": 2,
    "SO4": 2
}

am = Pitzer(ph)
print(am.get_lny_m("Mg"), 2.7182 ** am.get_lny_m("Mg"))
print(am.get_A())
print(am.get_harvie_j(1))
