from math import log, sqrt
import re
from pytherm import constants
from pytherm.activity.activity_model import Activity_model

R = constants.R


class PSRK:

    def __init__(self,
                 system: dict,
                 ms_params: dict,
                 activity_model: Activity_model) -> None:
        self.ms_params = ms_params
        self.system = system
        self.activity_model = activity_model

    def get_p(self, system, T, v):
        alphas = self.get_alphas(T)
        ai = self.get_ai(alphas)
        bi = self.get_bi()
        ge = self.get_ge(system, T)

        b = self.get_b(system, bi)
        a = self.get_a(system, T, b, ai, bi, ge)
        return (R * T) / (v - b) - a / v / (v + b)

    def get_alphas(self, T: float) -> float:
        # sqrt_T = sqrt(T)
        alphas = {}
        for i in self.system:
            Tr = T / self.ms_params[i]['Tc']
            sqrt_T = sqrt(Tr)
            if Tr < 1:
                c1 = self.ms_params[i]['c1']
                c2 = self.ms_params[i]['c2']
                c3 = self.ms_params[i]['c3']
                alphas[i] = (1
                             + c1 * (1 - sqrt_T)
                             + c2 * (1 - sqrt_T) ** 2
                             + c3 * (1 - sqrt_T) ** 3) ** 2
            else:
                c1 = self.ms_params[i]['c1']
                alphas[i] = (1
                             + c1 * (1 - sqrt_T)) ** 2
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

    def get_a(self, system, T, b, ai, bi, ge, A=-0.64663):
        s1 = 0
        s2 = 0
        for i in self.system:
            s1 += system[i] * ai[i] / bi[i]
            s2 += system[i] * log(b / bi[i])

        return b * (ge / A + s1 + R * T / A * s2)

    def get_b(self, system, bi):
        b = 0
        for i in self.system:
            b += system[i] * bi[i]
        return b

    def get_ge(self, system, T):
        if len(system) == 1:
            return 0
        else:
            return self.activity_model.get_ge(system, T)
