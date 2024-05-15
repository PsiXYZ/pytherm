from pytherm.systems.chemicalsystem import (
    ChemicalSystem,
    ChemicalReaction,
)
from pytherm.systems.phase import Phase, Solution, Substance
from pytherm.activity.electrolytemodel import ElectrolyteModel
import numpy as np


class EquilibriumSystem(ChemicalSystem):
    phases_list: list[Phase]
    reactions_list: list[ChemicalReaction]
    substances_list: list[Substance]
    T: float
    f_list: list[float]

    def __init__(self) -> None:
        self.phases_list = []
        self.reactions_list = []
        self.substances_list = []

    def add_phase(self, phase: Phase):
        if phase not in self.phases_list:
            self.phases_list.append(phase)

        for reaction in phase.reactions_list:
            if reaction not in self.reactions_list:
                self.reactions_list.append(reaction)

        for substance in phase.substances_list:
            if substance not in self.substances_list:
                self.substances_list.append(substance)

        self.f_list = np.full(len(self.reactions_list), 0.0)

    def equilibrate(
        self,
        T,
        k_lim=10,
        fabs=1e-5,
        ftol=1e-12,
    ):
        self.k_lim = k_lim
        self.fabs = fabs
        self.ftol = ftol

        self.T = T

        while 1:
            self.update_f()
            if np.sum(self.f_list) < self.fabs:
                break

            print("F = ", self.f_list)

            # поиск реакции с мах f
            max_i = 0
            max_f = -1
            for i, reaction in enumerate(self.reactions_list):
                if reaction.f > max_f:
                    max_i = i
                    max_f = reaction.f
            current_reaction = self.reactions_list[max_i]

            left_bound, right_bound = current_reaction.get_bounds()
            rs = np.array((left_bound, 0, right_bound))

            # начало оптимизации
            vals_u = np.array((0.0, 0.0))
            while 1:
                rs[1] = (rs[0] + rs[2]) / 2
                current_reaction.set_ksi(rs[1])

                current_reaction.phase.update_f()

                res = current_reaction.log10_k - current_reaction.log10_pr
                vals_u[1] = vals_u[0]
                vals_u[0] = current_reaction.f

                if res < 0:
                    rs[2] = rs[1]
                else:
                    rs[0] = rs[1]

                if vals_u[0] == 0:
                    current_reaction.set_ksi(rs[1])
                    current_reaction.fin_ksi()
                    break
                if vals_u[0] < self.fabs / 10:
                    current_reaction.set_ksi(rs[1])
                    current_reaction.fin_ksi()
                    break

                # tol = np.abs((vals_u[0] - vals_u[1]) / vals_u[0])
                # if tol < self.ftol:
                #     self.reactions_list[max_i].ksi = rs[1]
                #     break

    def update_f(self):
        for phase in self.phases_list:
            phase.update_f()

        for i, reaction in enumerate(self.reactions_list):
            self.f_list[i] = reaction.f
