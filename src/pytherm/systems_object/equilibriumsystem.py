from pytherm.systems_test.chemicalsystem import (
    ChemicalSystem,
    ChemicalReaction,
    Phase,
    Solution,
    SolidPhase,
)
from pytherm.activity.electrolytemodel import ElectrolyteModel
from pytherm.systems_test.equilibriumsolver import EquilibriumSolver
import numpy as np


class ElectrolyteSolution(Solution):
    activity_model: ElectrolyteModel

    def __init__(
        self,
        substances_dict: dict[str, float],
        reactions_list: list[ChemicalReaction],
        activity_model: ElectrolyteModel = None,
    ) -> None:
        super().__init__(
            substances_list=list(substances_dict.keys()),
            concentrations_list=list(substances_dict.values()),
            reactions_list=reactions_list,
        )
        self.activity_model = activity_model

    def set_activity_model(self, activity_model: ElectrolyteModel):
        self.activity_model = activity_model


class ExchangeInterface:
    pass


class SolubilitySP(ExchangeInterface):
    solid_phase: SolidPhase
    solution: Solution
    c0: np.ndarray
    log_k: np.ndarray

    def __init__(self, solid_phase: SolidPhase, solution: Solution) -> None:
        super().__init__()
        self.solid_phase = solid_phase
        self.solution = solution

    def set_reactions(self, reactions: list[ChemicalReaction]):
        pass


class EquilibriumSystem(ChemicalSystem):
    solver: EquilibriumSolver
    interfaces_list: list[ExchangeInterface]
    T: float

    def __init__(self) -> None:
        super().__init__()
        self.phases_list = []

    def add_phase(self, phase: Phase):
        self.phases_list.append(phase)

    def set_solver(self, solver: EquilibriumSolver):
        self.solver = solver

    def equilibrate(self, T):
        self.T = T
        for phase in self.phases_list:
            self.solver.equilibrate(phase)
