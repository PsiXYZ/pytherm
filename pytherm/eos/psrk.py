from math import log, sqrt
from pytherm import constants

R = constants.R


class PSRK:

    def __init__(self, system: dict = None, ms_params: dict = None) -> None:
        self.ms_params = ms_params
        self.system = system

    def get_p(self, system, T, v):
        alphas = self.get_alphas(T)
        ai = self.get_ai(alphas)
        bi = self.get_bi()
        ge = self.get_ge()

        b = self.get_b(system, bi)
        a = self.get_a(T, b, ai, bi, ge)
        b = 0
        return (R * T) / (v - b) - a / v / (v + b)

    def get_alphas(self, T: float) -> float:
        sqrt_T = sqrt(T)
        alphas = {}
        for i in self.system:
            c1 = self.ms_params[i]['c1']
            c2 = self.ms_params[i]['c2']
            c3 = self.ms_params[i]['c3']
            alphas[i] = (1 /
                         + c1 * (1 - sqrt_T) /
                         + c2 * (1 - sqrt_T) ** 2 /
                         + c3 * (1 - sqrt_T) ** 3) ** 2
        return alphas

    def get_ai(self, alphas):
        ai = {}
        for i in self.system:
            Tc = self.ms_params[i]['Tc']
            Pc = self.ms_params[i]['Pc']
            ai[i] = 0.42748 * R ** 2 * Tc ** 2 * alphas[i] / Pc
        return ai

    def get_bi(self):
        bi = {}
        for i in self.system:
            Tc = self.ms_params[i]['Tc']
            Pc = self.ms_params[i]['Pc']
            bi[i] = 0.08664 * R * Tc / Pc
        return bi

    def get_a(self, T, b, ai, bi, ge, A=-0.64663):
        s1 = 0
        s2 = 0
        for i in self.system:
            s1 += self.system[i] * ai[i] / bi[i]
            s2 += self.system[i] * log(b / bi[i])

        return b * (ge / A + s1 + R * T / A * s2)

    def get_b(self, system, bi):
        b = 0
        for i in self.system:
            b += self.system[i] * bi[i]
        return b

    def get_ge(self):
        return 0
