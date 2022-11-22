"""Module contain SIT class for activity calculations
"""
import math
from .electrolytes import get_I, get_A
from pytherm.stoichiometry import extract_charges
import numpy as np
from .db.sit.sit import sit_parameters

class SIT:
    substances = None
    charges = None
    cations = None
    anions = None
    epsilon_matrix = None

    def __init__(self, ph, parameters=sit_parameters):
        self.substances = np.array(list(ph.keys()))
        self.charges = np.array(extract_charges(ph))
        self.cations = self.substances[self.charges > 0]
        self.anions = self.substances[self.charges < 0]
        self.epsilon_matrix = np.zeros((len(self.cations), len(self.anions)))

        for i in range(len(self.cations)):
            for j in range(len(self.anions)):
                for k in parameters:
                    if k[0] == self.cations[i] and k[1] == self.anions[j]:
                        self.epsilon_matrix[i, j] = k[2]
        if .0 in self.epsilon_matrix:
            print("SIT model: zero values in epsilon_matrix")

    def get_y(self, ph: dict, T=298):
        molalities = np.array(list(ph.values()))

        A = get_A(T)
        I = get_I(molalities, charges=self.charges)
        sqrt_I = math.sqrt(I)
        D = A * sqrt_I / (1 + 1.5 * sqrt_I)

        lny = np.zeros((len(self.substances)))
        for i in range(len(self.substances)):
            if self.charges[i] > 0:
                j = np.where(self.cations == self.substances[i])[0][0]
                lny[i] = - self.charges[i] ** 2 * D + self.epsilon_matrix[j, :] @ molalities[self.charges < 0]
            else:
                j = np.where(self.anions == self.substances[i])[0][0]
                lny[i] = - self.charges[i] ** 2 * D + self.epsilon_matrix[:, j] @ molalities[self.charges > 0]

        y = {}
        for i in range(len(self.substances)):
            y[self.substances[i]] = np.exp(lny[i])
        return y
