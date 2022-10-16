from math import exp, log, sqrt
from pytherm import constants
from pytherm.activity.activity_model import Activity_model
import numpy as np
from .eos import EOS
R = constants.R


class PSRK(EOS):
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
        return (R * T) / (v - b) - a / (v * (v + b))

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
        return self.activity_model.get_ge(system, T)

    def get_cubic_coef(self, system, T, p):
        alphas = self.get_alphas(T)
        ai = self.get_ai(alphas)
        bi = self.get_bi()
        ge = self.get_ge(system, T)

        b = self.get_b(system, bi)
        a = self.get_a(system, T, b, ai, bi, ge)

        return (1,
                - R * T / p,
                a / p - b * R * T / p - b ** 2,
                - a * b / p)

    def get_roots(self, system, P, T):
        cs = self.get_cubic_coef(system=system, T=T, p=P)
        r = np.roots(cs)
        out = []
        for i in r:
            if i.imag == 0:
                out.append(i)
        if len(out) == 1:
            return [float(out[0])]
        else:
            return [float(min(r)), float(max(r))]

    def get_f(self, system, P, V, T):
        alphas = self.get_alphas(T)
        ai = self.get_ai(alphas)
        bi = self.get_bi()

        alp = {}
        for i in system:
            alp[i] = ai[i] / (bi[i] * R * T)

        b = self.get_b(system, bi)

        yi = self.activity_model.get_y(system, T)
        # s = 0
        # for i in system:
        #     s += system[i] * log(yi[i])
        ge_RT = self.activity_model.get_ge_RT(system, T)
        al = self.get_al(system, ge_RT, b, bi, alp)

        # a = self.get_a(system, T, b, ai, bi, ge)
        # roots = self.get_roots(system=system, P=P, T=T)
        der_ai = self.get_der_ai(yi, b, bi, alp)

        def f(v):
            lnf = {}
            for i in system:
                # lnf[i] = bi[i] / b * ((P * v) / (R * T) - 1)
                # - log(P * (v - b) / (R * T))
                # - der_ai[i] * log((v + b) / v)
                lnf[i] = log(R * T / (P * (v - b))) + (1 / (v - b) - al/(v + b)) * bi[i] - der_ai[i] * log((v + b) / v)
       
            fi = {}
            for i in lnf:
                fi[i] = exp(lnf[i])
            return fi
        return f(V)

    def get_der_ai(self, y, b, bi, alp, A=-0.64663):
        der_ai = {}
        for i in y:
            der_ai[i] = 1 / A * (log(y[i]) + log(b / bi[i]) +
                                 bi[i] / b - 1) + alp[i]
        return der_ai

    def get_al(self, system, ge_RT, b, bi, alp, A=-0.64663):
        s1 = 0
        s2 = 0
        for i in bi:
            s1 += system[i] * log(b / bi[i])
            s2 += system[i] * alp[i]
        return 1 / A * (ge_RT + s1) + s2

    def get_system(self):
        return self.system
