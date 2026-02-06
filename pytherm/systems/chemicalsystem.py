from __future__ import annotations

import numpy as np
from pytherm.stoichiometry import str_to_reaction
import pytherm.systems.phase as ph


class ChemicalReaction:
    log10_k: float
    reaction_rate: float
    reaction_vector: np.array
    substances_str_list: list[str]
    substances_list: list[ph.Substance]
    reaction_str: str
    phase: ph.Phase
    f: float
    ksi: float
    log10_pr: float

    def __init__(
        self,
        log10_k=None,
        reaction_vector=None,
        substances_str_list=None,
        reaction_str=None,
    ) -> None:
        self.log10_k = log10_k
        if reaction_str is None:
            self.reaction_vector = reaction_vector
            self.substances_str_list = substances_str_list
            self.reaction_str = reaction_str
        else:
            self.reaction_str = reaction_str
            self.reaction_vector, self.substances_str_list = str_to_reaction(
                reaction_str
            )

    def __str__(self) -> str:
        return self.reaction_str

    def get_bounds(self):      
        vals = []
        for i in range(len(self.reaction_vector)):
            if self.reaction_vector[i] > 0:
                buf = self.substances_list[i].c / self.reaction_vector[i]
                vals.append(-buf)
        r_l = max(vals)

        vals = []
        for i in range(len(self.reaction_vector)):
            if self.reaction_vector[i] < 0:
                buf = self.substances_list[i].c / self.reaction_vector[i]
                vals.append(-buf)
        r_r = min(vals)
        return r_l, r_r
    
    def set_ksi(self,  ksi):
        self.ksi = ksi
        for i, substance in enumerate(self.substances_list):
            substance.c = substance.c0 + self.reaction_vector[i] * ksi
        # self.phase.update_c()
    
    def fin_ksi(self):
        for i, substance in enumerate(self.substances_list):
            substance.c0 = substance.c
        self.ksi = 0
    
class ChemicalSystem:
    pass