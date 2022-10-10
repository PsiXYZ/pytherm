from math import log, sqrt


class PSRK:

    def __init__(self) -> None:
        self.ms_params = {}
        self.system = {}
        self.R = 8.314

    def get_p(self, T, P, system):
        pass

    def get_alphas(self, T: float) -> float:
        sqrt_T = sqrt(T)
        alphas = {}
        for i in self.system:
            c1 = self.ms_params[i]['c1']
            c2 = self.ms_params[i]['c2']
            c3 = self.ms_params[i]['c3']
            alphas[i] = (1 /
                         + c1[0] * (1 - sqrt_T) /
                         + c2[1] * (1 - sqrt_T) ** 2 /
                         + c3[2] * (1 - sqrt_T) ** 3) ** 2
        return alphas

    def get_ai(self, alphas):
        ai = {}
        for i in self.system:
            Tc = self.ms_params[i]['Tc']
            Pc = self.ms_params[i]['Pc']
            ai[i] = 0.42748 * self.R ** 2 * Tc ** 2 * alphas[i] / Pc
        return ai

    def get_bi(self, T, P):
        bi = {}
        for i in self.system:
            Tc = self.ms_params[i]['Tc']
            Pc = self.ms_params[i]['Pc']
            bi[i] = 0.08664 * self.R * Tc / Pc
        return bi

    def get_a(self, T, ai, bi, ge, A=-0.64663):
        b = 0
        for i in self.system:
            b += self.system[i] * bi[i]

        s1 = 0
        s2 = 0
        for i in self.system:
            s1 += self.system[i] * ai[i] / bi[i]
            s2 += self.system[i] * log(b / bi[i])

        return b * (ge / A + s1 + self.R * T / A * s2)

    def get_ge(self):
        return 1
