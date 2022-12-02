"""_summary_
"""

import numpy as np


class ChemicalReaction:
    def __init__(self,
                 log_k=None,
                 reaction_vector=None,
                 substances=None,
                 reaction_str=None
                 ) -> None:
        self.log_k = log_k
        self.reaction_vector = reaction_vector
        self.substances = substances
        self.reaction_str = reaction_str

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
            log_k.append(reaction.log_k)
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

    def find_equilibrium(self, ph):
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
        res = np.abs(np.log(pr) - self.log_k)
        return res

    def get_pr(self, ksi):
        pr = np.full(len(self.log_k), 0.0)  # реакционные произведения
        concentrations = self.get_concentrations(ksi)  # кол-ва вещества i

        for i in range(len(self.log_k)):
            activities = self.activity_model.get_a(concentrations)
            s = np.power(activities, self.reaction_matrix[i])
            pr[i] = np.prod(s)
        return pr

    def get_concentrations(self, ksi):
        concentrations = np.full(len(self.c_0), 0.0)
        for i in range(len(self.c_0)):
            concentrations[i] = self.reaction_matrix[:, i] @ ksi
        return self.c_0 + concentrations

    def get_bounds(self, ksi, reaction_index):
        p_ksi = ksi
        p_ksi[reaction_index] = 0
        c = self.get_concentrations(p_ksi)
        vals = []
        for i in range(len(self.c_0)):
            if self.reaction_matrix[reaction_index][i] > 0:
                buf = c[i] / self.reaction_matrix[reaction_index][i]
                vals.append(- buf)
        r_l = max(vals)

        vals = []
        for i in range(len(self.c_0)):
            if self.reaction_matrix[reaction_index][i] < 0:
                buf = c[i] / self.reaction_matrix[reaction_index][i]
                vals.append(- buf)
        r_r = min(vals)
        return r_l, r_r

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
