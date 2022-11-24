import numpy as np
from .db import pitzer as datasets
from pytherm.activity import electrolytes as el
from pytherm.activity.electrolytes import get_I
from pytherm.stoichiometry import extract_charges


class Pitzer:
    get_A: callable
    substances: np.ndarray
    charges: np.ndarray
    charges_dict: dict
    cations: np.ndarray
    anions: np.ndarray
    neutral: np.ndarray
    db: datasets.ParametersPitzerNew
    T_model = -1

    def __init__(self, ph, db=datasets.pitzer_dataset, get_A=el.get_A) -> None:
        self.get_A = get_A
        self.charges = np.array(extract_charges(list(ph.keys())))
        self.db = db

        subs = np.array(list(ph.keys()))
        self.cations = subs[self.charges > 0]
        self.anions = subs[self.charges < 0]
        self.neutral = subs[self.charges == 0]
        self.substances = subs
        self.charges_dict = {self.substances[i]: self.charges[i] for i in range(len(self.substances))}

    def get_y(self, ph: dict):
        y = {}
        for s in ph:
            y[s] = np.exp(self.grad(ph, s))
        return y

    def grad(self, ph, s, dm=1e-8):
        ph_b = {}
        for i in ph:
            ph_b[i] = ph[i]

        ph_b[s] = ph[s] - dm
        y1 = self.get_G_nRT(ph_b)
        ph_b[s] = ph[s] + dm
        y2 = self.get_G_nRT(ph_b)

        return (y2 - y1) / (2 * dm)

    def get_G_nRT(self, ph:dict, T=298.15):
        r"""Calculate :math:`G^{ex}/(w_{w}RT)`

        .. math::

            G^{ex}/(w_{w}RT) = f(I)
            + 2 \sum_{c}\sum_{a}[B_{ca}
            + (\sum_{c}m_c z_c)C_{ca}] \\
            + \sum_{c}\sum_{c'} m_c m_{c'} [2\Phi_{cc'} + \sum_a m_a \psi_{cc'a}] 
            + \sum_{a}\sum_{a'} m_a m_{a'} [2\Phi_{aa'} + \sum_c m_c \psi_{caa'}] \\
            + 2 \sum_n\sum_c m_n m_c \lambda_{nc}
            + 2 \sum_n\sum_a m_n m_a \lambda_{na} 
            + 2 \sum_n\sum_{n'} m_n m_{n'} \lambda_{nn'} + ... 

        Parameters
        ----------
        ph : _type_
            _description_
        """
        def get_ca():
            value_ca = 0
            for c in self.cations:
                m_c = ph[c]
                for a in self.anions:
                    if self.db.is_exist('CA', c, a):
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
                    if self.db.is_exist('THETA', c1, c2):
                        m_c2 = ph[c2]
                        fi = self.db['THETA'][c1][c2] \
                            + self.get_theta_e(c1, c2, I, A)
                        value_cc += m_c1 * m_c2 * 2 * fi
            return value_cc

        def get_aa():
            value_aa = 0
            for a1 in self.cations:
                m_a1 = ph[a1]
                for a2 in self.cations:
                    if self.db.is_exist('THETA', a1, a2):
                        m_a2 = ph[a2]
                        fi = self.db['THETA'][a1][a2]  \
                            + self.get_theta_e(a1, a2, I, A)
                        value_aa += m_a1 * m_a2 * 2 * fi
            return value_aa

        def get_cca():
            value_cca = 0
            for c1 in self.cations:
                m_c1 = ph[c1]
                for c2 in self.cations:
                    m_c2 = ph[c2]
                    for a in self.anions:
                        if self.db.is_exist('PSI', 'c1', 'c2', 'a'):
                            m_a = ph[a]
                            psi = self.db['PSI'][c1][c2][a]
                            value_cca += m_c1 * m_c2 * m_a * psi
            return value_cca

        def get_caa():
            value_caa = 0
            for a1 in self.cations:
                m_a1 = ph[a1]
                for a2 in self.cations:
                    m_a2 = ph[a2]
                    for c in self.anions:
                        if self.db.is_exist('PSI', c, a1, a2):
                            m_c = ph[c]
                            psi = self.db['PSI'][c][a1][a2]
                            value_caa += m_a1 * m_a2 * m_c * psi
            return value_caa

        def get_nc():
            value_nc = 0
            for n in self.neutral:
                m_n = ph[n]
                for c in self.cations:
                    m_c = ph[c]
                    if self.db.is_exist('LAMDA', n, c):
                        value_nc += 2 * m_n * m_c * self.db['LAMDA'][n][c]
            return value_nc

        def get_na():
            value_na = 0
            for n in self.neutral:
                m_n = ph[n]
                for a in self.anions:
                    m_a = ph[a]
                    if self.db.is_exist('LAMDA', n, a):
                        value_na += 2 * m_n * m_a * self.db['LAMDA'][n][a]
            return value_na

        def get_nca():
            value_nca = 0
            for n in self.cations:
                m_n = ph[n]
                for c in self.cations:
                    m_c = ph[c]
                    for a in self.anions:
                        if self.db.is_exist('ZETA', a, n, c):
                            m_a = ph[a]
                            value_nca += (m_n * m_c * m_a
                                          * self.db['ZETA'][n][c][a])
            return value_nca

        if self.T_model != T:
            self.db.update_params(T)
            self.T_model = T

        molalities = np.array(list(ph.values()))
        I = get_I(molalities, self.charges)
        A = self.get_A(T)
        Z = molalities @ np.abs(self.charges)

        G_nRT = self.get_f(I, A)
        G_nRT += get_ca()
        if len(self.cations) > 1:
            G_nRT += get_cc()
            G_nRT += get_cca()
        if len(self.anions) > 1:
            G_nRT += get_aa()
            G_nRT += get_caa()
        if len(self.neutral) > 0:
            G_nRT += get_nc()
            G_nRT += get_na()
            G_nRT += get_nca()

        return G_nRT

    def get_f(self, I, A, b=1.2):
        return - 4 * A * I * np.log(1 + b * np.sqrt(I)) / b

    def get_B(self, s1: str, s2: str, I: float):
        def get_g(x):
            return 2 * (1 - (1 + x) * np.exp(-x)) / x ** 2

        alpha = [0, 0]
        if np.abs(self.charges_dict[s1]) == 1 or np.abs(self.charges_dict[s2]) == 1:
            alpha = [2, 0]
        elif np.abs(self.charges_dict[s1]) == 2 and np.abs(self.charges_dict[s2]) == 2:
            alpha = [1.4, 12]
        else:
            alpha = [2, 50]

        b0 = self.db['CA'][s1][s2][0]
        b1 = self.db['CA'][s1][s2][1]
        b2 = self.db['CA'][s1][s2][2]
        g1 = get_g(alpha[0] * np.sqrt(I))
        if alpha[1] == 0:
            return b0 + b1 * g1
        else:
            g2 = get_g(alpha[1] * np.sqrt(I))
            return b0 + b1 * g1 + b2 * g2

    def get_C(self, s1: str, s2: str):
        C_f = self.db['CA'][s1][s2][3]
        z1 = self.charges_dict[s1]
        z2 = self.charges_dict[s2]
        return C_f / (2 * np.sqrt(np.abs(z1 * z2)))

    def get_theta_e(self, s1: str, s2: str, I: float, A: float):
        r"""Calculate excess term :math:`^E\theta_{ij}` for :math:`\Phi_{ij}` using

        .. math::
            \theta_{MN}=\frac{z_M z_N}{4I}[J(x_{MN})- 0.5 J(x_{MM} )- 0.5 J(x_{NN})]

        Parameters
        ----------
        s1 : _type_
            _description_
        s2 : _type_
            _description_
        I : _type_
            _description_

        Returns
        -------
        _type_
            _description_
        """
        def get_x_ij(z1, z2, I, A):
            return 6 * z1 * z2 * A * np.sqrt(I)

        z_m = self.charges_dict[s1]
        z_n = self.charges_dict[s2]

        x_mn = get_x_ij(z_m, z_n, I, A)
        x_mm = get_x_ij(z_m, z_m, I, A)
        x_nn = get_x_ij(z_n, z_n, I, A)
        J1, Jp1 = self.get_harvie_j(x_mn)
        J2, Jp2 = self.get_harvie_j(x_mm)
        J3, Jp3 = self.get_harvie_j(x_nn)

        E = (z_m * z_n / (4 * I)) * (J1 - 0.5 * J2 - 0.5 * J3)
        # Ep = (- E / I + (z_m * z_n / (8 * I ** 2))
        #       * (x_mn * Jp1 - 0.5 * x_mm * Jp2 - 0.5 * x_nn * Jp3))
        return E

    def get_harvie_j(self, x):
        """Calculate J using Chebyshev polynomial approximations

        Parameters
        ----------
        x : _type_
            _description_

        Returns
        -------
        _type_
            _description_
        """
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