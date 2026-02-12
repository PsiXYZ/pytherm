from __future__ import annotations

import numpy as np
from pytherm.stoichiometry import str_to_reaction
import re
import scipy


class EquilibriumReaction:
    """Class for chemical reaction representation"""

    log10_k: float
    reaction_rate: float
    reaction_vector: np.array
    substances_list: list[str]
    reaction_str: str

    def __init__(
        self,
        log10_k=None,
        reaction_vector=None,
        substances_list=None,
        reaction_str=None,
    ) -> None:
        self.log10_k = log10_k
        if reaction_str is None:
            self.reaction_vector = reaction_vector
            self.substances_list = substances_list
            self.reaction_str = reaction_str
        else:
            self.reaction_str = reaction_str
            self.reaction_vector, self.substances_list = str_to_reaction(reaction_str)

    def __str__(self) -> str:
        return self.reaction_str

class Phase:
    substances_list: list[str]

    def __init__(self, substances_list) -> None:
        self.substances_list = substances_list.keys()
        self.c = substances_list.values()

def make_reaction_matrix(reactions: list[EquilibriumReaction]):
    substances_list = []
    for reaction in reactions:
        for substance in reaction.substances_list:
            if substance not in substances_list:
                substances_list.append(substance)

    reaction_matrix = np.zeros((len(reactions), len(substances_list)))
    for reaction_i in range(len(reactions)):
        for substance_i in range(len(reactions[reaction_i].substances_list)):
            i = substances_list.index(
                reactions[reaction_i].substances_list[substance_i]
            )
            reaction_matrix[reaction_i, i] = reactions[reaction_i].reaction_vector[
                substance_i
            ]
    return substances_list, reaction_matrix


class EquilibriumSystem():
    T: float
    reactions_v: list[EquilibriumReaction]
    phases_v: list[Phase]
    substances_list: list[str]
    phases_names: list[str]
    substances_raw: list[str]
    
    log10_k: np.ndarray

    Ar: np.ndarray
    Af: np.ndarray
    n0: np.ndarray
    n: np.ndarray
    n_tot: np.ndarray
    c: np.ndarray
    a: np.ndarray
    f: np.ndarray

    def __init__(self) -> None:
        super().__init__()
        self.phases_v = []
        self.substances_list = []
        self.phases_names = []
        self.substances_raw = []

    def init_by_reactions(self, reactions_list: list[EquilibriumReaction]):
        self.reactions_v = reactions_list
        r = make_reaction_matrix(reactions_list)
        self.substances_raw = r[0]
        self.Ar = r[1]

        self.substances_list = []
        self.phases_names = []
        for substance in self.substances_raw:
            self.substances_list.append(substance.split("{")[0])
            if substance.split("{")[1][:-1] not in self.phases_names:
                self.phases_names.append(substance.split("{")[1][:-1])
        
        self.Af = np.zeros((len(self.phases_names), len(self.substances_list)))
        for substance in self.substances_raw:
            for phase in self.phases_names:
                if phase in substance:
                    self.Af[self.phases_names.index(phase), self.substances_list.index(substance.split("{")[0])] = 1




    def add_phase_by_dict(self, phase: dict[str, float]):
        names = list(phase.keys())
        phase_name = re.search(r"\{([^}]*)\}", names[0]).group(1)

        if phase_name not in self.phases_names:
            self.phases_names.append(phase_name)
        
        for substance in names:
            substance_name = re.sub(r"\{[^}]*\}", "", substance)
            self.substances_list.append(substance_name)
            self.substances_raw.append(substance)
              
    def add_reactions_by_list(self, reactions_list: list[EquilibriumReaction]):
        self.reactions_v = reactions_list

    def assemble(self):
        self.Ar = np.zeros((len(self.reactions_v), len(self.substances_list)))
        self.Af = np.zeros((len(self.phases_names), len(self.substances_list)))

        r = make_reaction_matrix(self.reactions_v)
        for i in range(len(self.substances_raw)):
            if self.substances_raw[i] in r[0]:
                self.Ar[:, i] = r[1][:, r[0].index(self.substances_raw[i])]
        
        for i in range(len(self.phases_names)):
            for j in range(len(self.substances_raw)):
                if self.phases_names[i] in self.substances_raw[j]:
                    self.Af[i, j] = 1
        
        self.log10_k = np.zeros(len(self.reactions_v))
        for i in range(len(self.reactions_v)):
            self.log10_k[i] = self.reactions_v[i].log10_k

        self.n0 = np.zeros(len(self.substances_list))
        self.n = np.zeros(len(self.substances_list))
        self.Ar = self.Ar.T

    def get_Q(self, ksi: np.ndarray) -> np.ndarray:
        self.n = self.n0 + self.Ar @ ksi
        self.n_tot = self.Af @ self.n 
        self.c = self.n / (self.Af.T @ self.n_tot)
        self.a = self.c

        # lg_a = np.log10(self.a)
        # lg_a = np.where(np.isinf(lg_a), 0.0, lg_a)
        # p = lg_a @ self.Ar
        # p = np.where(np.isnan(p), np.inf, p)

        a_safe = np.clip(self.a, 1e-300, None)
        Q = np.prod(a_safe[:, None] ** self.Ar, axis=0)
        # Q = np.prod(np.power(self.a[:, None], self.Ar), axis=0)
        # lg_Q = np.log10(Q)
        return Q
    
    # def get_f(self, ksi: np.ndarray) -> np.ndarray:
    #     f = self.get_Q(ksi) - self.log10_k
    #     # f = np.where(np.isnan(f), np.inf, f)
    #     return f

    def get_c(self, ksi: np.ndarray) -> np.ndarray:
        self.n = self.n0 + self.Ar @ ksi
        self.n_tot = self.Af @ self.n 
        self.c = self.n / (self.Af.T @ self.n_tot)
        return self.c
    
    def get_n(self, ksi: np.ndarray) -> np.ndarray:
        self.n = self.n0 + self.Ar @ ksi
        return self.n

    def get_ni(self, ksi: np.ndarray, i: int) -> float:
        n = self.get_n(ksi)
        return n[i]

    def get_bounds(self, ksi: np.ndarray, reaction_index: int):
        p_ksi = ksi.copy()
        p_ksi[reaction_index] = 0
        n = self.get_n(p_ksi)

        vals = - n / self.Ar[:, reaction_index]
        vals = np.where(np.isinf(vals), 0.0, vals)
        vals = np.where(np.isnan(vals), 0.0, vals)
        r_l = max(vals)
        r_r = min(vals)

        # vals = []
        # for i in range(len(self.Ar[0])):
        #     if self.Ar[reaction_index, i] > 0:
        #         buf = c[i] / self.Ar[reaction_index, i]
        #         vals.append(-buf)
        # r_l = max(vals)

        # vals = []
        # for i in range(len(self.Ar[0])):
        #     if self.Ar[reaction_index, i] < 0:
        #         buf = n[i] / self.Ar[reaction_index, i]
        #         vals.append(-buf)
        # r_r = min(vals)
        return r_l, r_r

    def equilibrate(self, ph, solver_type='scipy'):
        ph_flat = {}
        for i in ph:
            for j in ph[i]:
                ph_flat[j] = ph[i][j]
        
        for i in ph_flat:
            self.n0[self.substances_raw.index(i)] = ph_flat[i]

        if solver_type == 'scipy':
            solver = ScipyEquilibriumSolver()
        elif solver_type == 'custom_lm':
            solver = LevenbergMarquardtSolver()
        else:
            solver = EquilibriumSolver()

        self.ksi = solver.equilibrate(self)


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

    def equilibrate(self, system: EquilibriumSystem) -> np.ndarray:
        n_reactions = len(system.log10_k)
        ksi = np.full(n_reactions, 0, dtype=float)
        
        while (1):
            F = np.abs(np.log10(system.get_Q(ksi)) - system.log10_k)
            # print("F = ", F)
            # print("ksi = ", ksi)

            if np.sum(np.abs(F)) < self.fabs:
                self.ksi = ksi
                break

            # поиск реакции с мах f
            v = 0
            ri = 0
            for i in range(len(F)):
                if F[i] > v:
                    v = F[i]
                    ri = i

            p_ksi = ksi
            p_ksi[ri] = 0
            left_bound, right_bound = system.get_bounds(p_ksi, ri)
            rs = np.array((left_bound, 0, right_bound))

            # начало оптимизации
            vals_u = np.array((0.0, 0.0))
            while (1):
                rs[1] = (rs[0] + rs[2]) / 2
                p_ksi[ri] = rs[1]
                # print("F = ", fi(p_ksi))
                # print("Ksi = ", p_ksi)
                # pr = system.get_p(p_ksi)
                # log10_K = system.log10_k

                c = system.get_c(p_ksi)
                Q = system.get_Q(p_ksi)
                log10_Q = np.log10(Q)
                
                res = (system.log10_k - np.log10(system.get_Q(p_ksi)))
                vals_u[1] = vals_u[0]
                vals_u[0] = np.abs(system.get_Q(p_ksi)[ri] - system.log10_k[ri])

                if res[ri] > 0:
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
        
        return ksi

class ScipyEquilibriumSolver:
    def __init__(self, method='lm', tol=1e-12):
        self.method = method
        self.tol = tol

    def equilibrate(self, system: EquilibriumSystem) -> np.ndarray:
        n_reactions = len(system.log10_k)
        ksi0 = np.full(n_reactions, 0.0, dtype=float)

        def residuals(ksi):
            Q = system.get_Q(ksi)
            # Avoid log(0) or very small numbers
            Q_safe = np.where(Q <= 1e-300, 1e-300, Q)
            er = np.abs(np.log10(Q_safe) - system.log10_k)
            # print(er)
            return np.sum(er)
        
        def residuals2(ksi):
            Q = system.get_Q(ksi)
            # Avoid log(0) or very small numbers
            Q_safe = np.where(Q <= 1e-300, 1e-300, Q)
            er = np.abs(np.log10(Q_safe) - system.log10_k)
            # print(er)
            return er

        sol = scipy.optimize.root(residuals2, ksi0, method="lm", tol=1e-6)
        # sol = scipy.optimize.minimize(residuals, ksi0, method='Nelder-Mead', tol=1e-12, options={'maxiter': 10000})

        # constraints = []
        # for i in range(len(system.substances_list)):
        #     constraints.append({'type': 'ineq', 'fun': lambda ksi, i=i: system.get_ni(ksi, i)})
        
        # constraints = []
        # for i in range(len(system.substances_list)):
        #     constraints.append(scipy.optimize.NonlinearConstraint(lambda ksi, i=i: system.get_ni(ksi, i), 0, np.inf))

        # print(constraints[9].fun(ksi0))

        # sol = scipy.optimize.minimize(residuals, ksi0, method='cobyla', tol=1e-12, options={'maxiter': 10000}, 
        #     constraints=constraints,
        # )

        # bounds = []
        # for i in range(len(system.reactions_v)):
        #     bounds.append((0, 0.3))    
        # sol = scipy.optimize.direct(residuals, bounds, len_tol=1e-9, maxiter=10000)
        # # sol = scipy.optimize.shgo(residuals, bounds)
        

        
        # Ensure final state is set
        system.get_Q(sol.x)

        if not sol.success:
            print(f"Solver failed: {sol.message}")
        
        return sol.x


class LevenbergMarquardtSolver:
    def __init__(self, atol=1e-6, tol=1e-6, max_iter=100, lambda_init=0.01, lambda_factor=10.0):
        self.atol = atol
        self.tol = tol
        self.max_iter = max_iter
        self.lambda_init = lambda_init
        self.lambda_factor = lambda_factor

    def _jacobian(self, system, ksi, residuals_val):
        """
        Compute numerical Jacobian of residuals with respect to ksi using finite differences.
        residuals_val is the current value of residuals(ksi).
        """
        n = len(ksi)
        m = len(residuals_val)
        J = np.zeros((m, n))
        epsilon = 1e-8  # Small perturbation

        for i in range(n):
            ksi_perturbed = ksi.copy()
            ksi_perturbed[i] += epsilon
            
            # Recalculate residuals for perturbed ksi
            Q_perturbed = system.get_Q(ksi_perturbed)
            Q_safe_p = np.where(Q_perturbed <= 1e-300, 1e-300, Q_perturbed)
            residuals_p = np.log10(Q_safe_p) - system.log10_k
            
            J[:, i] = (residuals_p - residuals_val) / epsilon
        
        return J

    def equilibrate(self, system: EquilibriumSystem) -> np.ndarray:
        n_reactions = len(system.log10_k)
        ksi = np.full(n_reactions, 0.0, dtype=float)
        lambda_val = self.lambda_init

        for iteration in range(self.max_iter):
            # Calculate residuals
            Q = system.get_Q(ksi)
            Q_safe = np.where(Q <= 1e-300, 1e-300, Q)
            residuals = np.log10(Q_safe) - system.log10_k
            
            # Check for convergence
            cost = np.sum(residuals**2)
            # cost = np.sum(np.abs(residuals))
            # if cost < self.tol:
            #     # Converged
            #     break

            if np.sum(np.abs(residuals)) < self.atol:
                # Converged
                break

            # Calculate Jacobian
            J = self._jacobian(system, ksi, residuals)

            # JTJ + lambda * I * diag(JTJ)
            JTJ = J.T @ J
            
            # Standard Levenberg-Marquardt damping
            H = JTJ + lambda_val * np.eye(n_reactions)

            g = J.T @ residuals

            try:
                delta = np.linalg.solve(H, -g)
            except np.linalg.LinAlgError:
                # If singular, increase damping and try again (simple retry logic)
                lambda_val *= self.lambda_factor
                continue

            # Evaluate new point
            ksi_new = ksi + delta
            
            Q_new = system.get_Q(ksi_new)
            Q_safe_new = np.where(Q_new <= 1e-300, 1e-300, Q_new)
            residuals_new = np.log10(Q_safe_new) - system.log10_k
            cost_new = np.sum(residuals_new**2)

            if cost_new < cost:
                # Accept step
                ksi = ksi_new
                lambda_val /= self.lambda_factor
                
                # Check absolute convergence on change
                if np.linalg.norm(delta) < self.tol and (cost - cost_new) < self.tol:
                    #  break
                    pass
            else:
                # Reject step and increase damping
                lambda_val *= self.lambda_factor

        system.ksi = ksi
        return ksi
