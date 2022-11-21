"""_summary_
"""

import stoichiometry as sm
import numpy as np


class ChemicalReaction:
    def __init__(self,
                 K_eqilibrium=None,
                 reaction_vector=None,
                 substances=None,
                 reaction_str=None
                 ) -> None:
        self.K_eqilibrium = K_eqilibrium
        self.reaction_vector = reaction_vector
        self.substances = substances
        self.reaction_str = sm.reaction_to_str(reaction_vector, substances)

    def __str__(self) -> str:
        return self.reaction_str


class EqilibriumSystem:
    substances = None
    reactions: list[ChemicalReaction]
    reaction_matrix: np.ndarray

    insol = None
    solubilities: list[ChemicalReaction]
    solubilities_matrix: np.ndarray

    def set_reactions(self, reactions: list[ChemicalReaction]):
        self.reactions = reactions
        self.substances, self.reaction_matrix = make_reaction_matrix(
            reactions
        )

    def set_solubilities(self, reactions: list[ChemicalReaction]):
        self.solubilities = reactions
        self.insol, self.solubilities_matrix = make_reaction_matrix(
            reactions
        )


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


# Ca_+2{aq} Na_+1{aq} CaSO4{s}
reactions = [
    ChemicalReaction(
        reaction_vector=[-1, 1, 1],
        substances=['H2O', 'H', 'OH']
    ),
    ChemicalReaction(
        reaction_vector=[-1, -1, 1],
        substances=['Ca', 'SO4', 'CaSO4']
    ),
    ChemicalReaction(
        reaction_vector=[-1, 1, 1],
        substances=['HSO4', 'H', 'SO4']
    ),
]

s = EqilibriumSystem()
s.set_reactions(reactions=reactions)
