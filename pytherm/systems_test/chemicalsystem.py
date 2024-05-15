import numpy as np
from pytherm.stoichiometry import str_to_reaction
from pytherm.activity.activitymodel import ActivityModel


class ChemicalReaction:
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
    concentrations_list: list[float]
    reactions_list: list[ChemicalReaction]
    log_k: np.ndarray
    reaction_matrix: np.ndarray

    def __init__(self) -> None:
        pass

    def fi(self, ksi: np.ndarray):
        return None

    def get_pr(self, ksi: np.ndarray):
        return None

    def get_concentrations(self, ksi: np.ndarray):
        return None


class Solution(Phase):
    activity_model: ActivityModel
    volume: float
    density: float

    def __init__(
        self,
        substances_list: list[str],
        concentrations_list: np.ndarray,
        reactions_list: list[ChemicalReaction],
    ) -> None:
        super().__init__()
        self.volume = 1
        self.density = 1

        self.set_reactions(reactions_list)
        for substance in substances_list:
            if not substance in self.substances_list:
                self.substances_list.append(substance)

        self.concentrations_list = np.full(len(self.substances_list), 0.0)
        for i1, s1 in enumerate(substances_list):
            for i2, s2 in enumerate(self.substances_list):
                if s1 == s2:
                    self.concentrations_list[i2] = concentrations_list[i1]

    def set_reactions(self, reactions: list[ChemicalReaction]):
        self.reactions = reactions
        self.substances_list, self.reaction_matrix = make_reaction_matrix(reactions)

        log_k = []
        for reaction in self.reactions:
            log_k.append(reaction.log10_k)
        self.log_k = np.array(log_k)

    def fi(self, ksi: np.ndarray):
        pr = self.get_pr(ksi)
        res = np.abs(np.log10(pr) - self.log_k)
        return res

    def get_pr(self, ksi):
        pr = np.full(len(self.log_k), 0.0)  # реакционные произведения
        concentrations = self.get_concentrations(ksi)  # кол-ва вещества i
        activities = self.activity_model.get_a(concentrations)
        for i in range(len(self.log_k)):
            s = np.power(
                activities[: len(self.reaction_matrix[0])], self.reaction_matrix[i]
            )
            for j in range(len(s)):
                if s[j] == np.inf:
                    s[j] = 0
            pr[i] = np.prod(s)
        return pr

    def get_concentrations(self, ksi):
        concentrations = np.full(len(self.concentrations_list), 0.0)
        for i in range(len(self.reaction_matrix[0])):
            concentrations[i] = self.reaction_matrix[:, i] @ ksi
        return self.concentrations_list + concentrations

    def get_bounds(self, ksi, reaction_index):
        p_ksi = ksi
        p_ksi[reaction_index] = 0
        c = self.get_concentrations(p_ksi)
        vals = []
        for i in range(len(self.reaction_matrix[0])):
            if self.reaction_matrix[reaction_index][i] > 0:
                buf = c[i] / self.reaction_matrix[reaction_index][i]
                vals.append(-buf)
        r_l = max(vals)

        vals = []
        for i in range(len(self.reaction_matrix[0])):
            if self.reaction_matrix[reaction_index][i] < 0:
                buf = c[i] / self.reaction_matrix[reaction_index][i]
                vals.append(-buf)
        r_r = min(vals)
        return r_l, r_r


class SolidPhase(Phase):
    def __init__(self, substances_list) -> None:
        super().__init__()
        
        self.substances_list = substances_list


class GasPhase(Phase):
    pass


class ChemicalSystem:
    phases_list: list[Phase]
    substances_list: list[str]


def make_reaction_matrix(reactions: list[ChemicalReaction]):
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
