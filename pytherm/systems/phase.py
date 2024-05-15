from __future__ import annotations

import pytherm.systems.chemicalsystem as chemsys
from pytherm.activity.activitymodel import ActivityModel
import numpy as np


class Substance:
    name: str
    n: float
    c: float
    c0: float
    a: float

    def __init__(self, name) -> None:
        self.name = name
        self.a = 0
        self.c = 0
        self.c0 = 0
        self.n = 0


class Phase:
    substances_list: list[Substance]
    substances_str_list: list[str]
    reactions_list: list[chemsys.ChemicalReaction]

    def __init__(self) -> None:
        pass

    def update_f(self):
        pass

    def update_pr(self):
        pass


class Solution(Phase):
    activity_model: ActivityModel
    volume: float
    density: float
    f_list: np.ndarray
    conc: np.ndarray
    ksi: np.ndarray

    def __init__(self) -> None:
        super().__init__()
        self.substances_str_list = []
        self.substances_list = []
        self.volume = 1

    def add_reactions(self, reactions_list: list[chemsys.ChemicalReaction]):
        self.reactions_list = reactions_list

        for reaction in self.reactions_list:
            for substance in reaction.substances_str_list:
                if substance not in self.substances_str_list:
                    self.substances_str_list.append(substance)
                    self.substances_list.append(Substance(substance))

        for reaction in self.reactions_list:
            reaction.phase = self

            reaction.substances_list = []
            for substance1 in reaction.substances_str_list:
                for i, substance2 in enumerate(self.substances_str_list):
                    if substance1 == substance2:
                        reaction.substances_list.append(self.substances_list[i])

        self.conc = np.full(len(self.substances_list), 0.0)

    def add_substances(self, substances_str_list: list[str]):
        for substance in substances_str_list:
            if substance not in self.substances_str_list:
                self.substances_str_list.append(substance)
                self.substances_list.append(Substance(substance))

        self.conc = np.full(len(self.substances_list), 0.0)

    def set_activity_model(self, activity_model: ActivityModel):
        self.activity_model = activity_model

    def update_f(self):
        self.update_pr()
        for reaction in self.reactions_list:
            reaction.f = abs(reaction.log10_pr - reaction.log10_k)

    def update_a(self):
        for i, substance in enumerate(self.substances_list):
            self.conc[i] = substance.c

        a = self.activity_model.get_a(self.conc)
        for i, substance in enumerate(self.substances_list):
            self.substances_list[i].a = a[i]

    def update_pr(self):
        self.update_a()

        for reaction in self.reactions_list:
            pr = 1
            for i, substance in enumerate(reaction.substances_list):
                pr *= substance.a ** reaction.reaction_vector[i]
            reaction.log10_pr = np.log10(pr)


class ElectrolyteSolution(Solution):
    def __init__(self) -> None:
        super().__init__()

    def set_concentrations(self, ph: dict[str, float]):
        for i, substance in enumerate(self.substances_str_list):
            if substance in ph:
                self.substances_list[i].c0 = ph[substance]
                self.substances_list[i].c = ph[substance]
