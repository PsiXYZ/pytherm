"""This module contains classes and functions for equilibrium calculations
"""

import numpy as np
from .stoichiometry import str_to_reaction


class ChemicalReaction:
    """Class for chemical reaction representation
    """
    def __init__(self,
                 log10_k=None,
                 reaction_vector=None,
                 substances=None,
                 reaction_str=None
                 ) -> None:

        self.log10_k = log10_k
        if reaction_str is None:
            self.reaction_vector = reaction_vector
            self.substances = substances
            self.reaction_str = reaction_str
        else:
            self.reaction_str = reaction_str
            self.reaction_vector, self.substances = str_to_reaction(reaction_str)

    def __str__(self) -> str:
        return self.reaction_str


class EquilibriumSystem:
    substances = None
    reactions: list[ChemicalReaction]
    reaction_matrix: np.ndarray
    log_k = None

    insol = None
    solids: list[ChemicalReaction]
    solids_matrix: np.ndarray

    activity_model = None

    c_0 = None

    def __init__(self, subs, reactions):
        self.set_reactions(reactions)
        for s in subs:
            if not s in self.substances:
                self.substances.append(s)

    def set_reactions(self, reactions: list[ChemicalReaction]):
        self.reactions = reactions
        substances, self.reaction_matrix = make_reaction_matrix(
            reactions
        )
        if self.substances is not None:
            subs = self.substances
            for s in subs:
                if s not in substances:
                    substances.append(s)
            self.substances = substances
        else:
            self.substances = substances

        log_k = []
        for reaction in self.reactions:
            log_k.append(reaction.log10_k)
        self.log_k = np.array(log_k)

    def set_solids(self, reactions: list[ChemicalReaction]):
        self.solids = reactions
        self.insol, self.solids_matrix = make_reaction_matrix(
            reactions
        )

    def set_solutes(self, solutes):
        self.substances = solutes

    def set_activity_model(self, activity_model):
        self.activity_model = activity_model

    def set_concentrations(self, ph):
        concentrations_dict = {}
        for s in self.substances:
            concentrations_dict[s] = 0
        for s in ph:
            concentrations_dict[s] = ph[s]
        concentrations = np.array(list(concentrations_dict.values()))
        self.c_0 = concentrations

        print(concentrations)

    def fi(self, ksi):
        pr = self.get_pr(ksi)
        res = np.abs(np.log10(pr) - self.log_k)
        return res

    def get_pr(self, ksi):
        pr = np.full(len(self.log_k), 0.0)  # реакционные произведения
        concentrations = self.get_concentrations(ksi)  # кол-ва вещества i
        # conc_dict = self.conc_to_dict(concentrations)
        activities = self.activity_model.get_a(concentrations)
        for i in range(len(self.log_k)):

            s = np.power(activities[:len(self.reaction_matrix[0])], self.reaction_matrix[i])
            for j in range(len(s)):
                if s[j] == np.inf:
                    s[j] = 0
            pr[i] = np.prod(s)
        return pr

    def get_concentrations(self, ksi):
        concentrations = np.full(len(self.c_0), 0.0)
        for i in range(len(self.reaction_matrix[0])):
            concentrations[i] = self.reaction_matrix[:, i] @ ksi
        # for i in range(len(self.reaction_matrix), len(self.c_0)):
        #     concentrations[i] = 0
        return self.c_0 + concentrations

    def get_bounds(self, ksi, reaction_index):
        p_ksi = ksi
        p_ksi[reaction_index] = 0
        c = self.get_concentrations(p_ksi)
        vals = []
        for i in range(len(self.reaction_matrix[0])):
            if self.reaction_matrix[reaction_index][i] > 0:
                buf = c[i] / self.reaction_matrix[reaction_index][i]
                vals.append(- buf)
        r_l = max(vals)

        vals = []
        for i in range(len(self.reaction_matrix[0])):
            if self.reaction_matrix[reaction_index][i] < 0:
                buf = c[i] / self.reaction_matrix[reaction_index][i]
                vals.append(- buf)
        r_r = min(vals)
        return r_l, r_r

    def conc_to_dict(self, conc):
        d = {}
        for i in range(len(self.substances)):
            d[self.substances[i]] = conc[i]
        return d


class EquilibriumSolver:
    eq_sys: EquilibriumSystem
    k_lim: float
    fabs: float
    ftol: float
    ksi: list[float]

    def __init__(self,
                 eq_sys: EquilibriumSystem,
                 k_lim=10,
                 fabs=1e-5,
                 ftol=1e-12,
                 ):
        self.eq_sys = eq_sys
        self.k_lim = k_lim
        self.fabs = fabs
        self.ftol = ftol

    def equilibrate(self):
        n_reactions = len(self.eq_sys.log_k)
        for i in range(n_reactions):
            pass

        ksi = np.full(n_reactions, 0, dtype=float)
        while (1):
            f = self.eq_sys.fi(ksi)
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
            left_bound, right_bound = self.eq_sys.get_bounds(p_ksi, ri)
            rs = np.array((left_bound, 0, right_bound))

            # начало оптимизации
            vals_u = np.array((0.0, 0.0))
            while (1):
                rs[1] = (rs[0] + rs[2]) / 2
                p_ksi[ri] = rs[1]
                # print("F = ", fi(p_ksi))
                # print("Ksi = ", p_ksi)
                pr = self.eq_sys.get_pr(p_ksi)
                K = self.eq_sys.log_k

                res = (K - np.log10(pr))
                vals_u[1] = vals_u[0]
                vals_u[0] = self.eq_sys.fi(p_ksi)[ri]

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


def make_reaction_matrix(reactions: list[ChemicalReaction]):
    substances = []
    for reaction in reactions:
        for substance in reaction.substances:
            if substance not in substances:
                substances.append(substance)

    reaction_matrix = np.zeros((len(reactions), len(substances)))
    for reaction_i in range(len(reactions)):
        for substance_i in range(len(reactions[reaction_i].substances)):
            i = substances.index(
                reactions[reaction_i].substances[substance_i])
            reaction_matrix[reaction_i,
                            i] = reactions[reaction_i].reaction_vector[substance_i]
    return substances, reaction_matrix
