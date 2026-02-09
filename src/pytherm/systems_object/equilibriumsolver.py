from pytherm.systems_test.equilibriumsystem import Phase
import numpy as np


class EquilibriumSolver:
    k_lim: float
    fabs: float
    ftol: float
    ksi: list[float]

    def __init__(
        self,
        k_lim=10,
        fabs=1e-5,
        ftol=1e-12,
    ):
        self.k_lim = k_lim
        self.fabs = fabs
        self.ftol = ftol

    def equilibrate(self, phase: Phase):
        n_reactions = len(phase.log_k)
        ksi = np.full(n_reactions, 0, dtype=float)
        
        while (1):
            f = phase.fi(ksi)
            print("F = ", f)
            print("ksi = ", ksi)

            if np.sum(f) < self.fabs:
                self.ksi = ksi
                break

            # поиск реакции с мах f
            v = 0
            ri = 0
            for i in range(len(f)):
                if f[i] > v:
                    v = f[i]
                    ri = i

            p_ksi = ksi
            p_ksi[ri] = 0
            left_bound, right_bound = phase.get_bounds(p_ksi, ri)
            rs = np.array((left_bound, 0, right_bound))

            # начало оптимизации
            vals_u = np.array((0.0, 0.0))
            while (1):
                rs[1] = (rs[0] + rs[2]) / 2
                p_ksi[ri] = rs[1]
                # print("F = ", fi(p_ksi))
                # print("Ksi = ", p_ksi)
                pr = phase.get_pr(p_ksi)
                K = phase.log_k

                res = (K - np.log10(pr))
                vals_u[1] = vals_u[0]
                vals_u[0] = phase.fi(p_ksi)[ri]

                if res[ri] < 0:
                    rs[2] = rs[1]
                else:
                    rs[0] = rs[1]

                if vals_u[0] == 0:
                    ksi[ri] = rs[1]
                    break
                if vals_u[0] < self.fabs / 10:
                    ksi[ri] = rs[1]
                    break
                tol = np.abs((vals_u[0] - vals_u[1]) / vals_u[0])
                if tol < self.ftol:
                    ksi[ri] = rs[1]
                    break