"""Module contain SIT class for activity calculations

.. contents:: :local:

SIT Class
------------
.. autoclass:: SIT
    :members: get_y
    :undoc-members:
    :show-inheritance:
    :member-order: bysource
"""
import math
from .electrolytes import get_I, get_A
from pytherm.stoichiometry import extract_charges
import numpy as np
from .db.sit.sit import sit_parameters


class SIT:
    r""" Class for activity calculations using Specific ion interaction theory

    .. math::
        \log_{10}(\gamma_{i}) = - z_i^2 * \frac{0.51 * \sqrt{I}}{1 + 1.5\sqrt{I}} + \sum_k\epsilon_{ik}m_k

    Parameters
    ----------
    ph: dict
        Phase dict for parameters initialization
    parameters : list
        SIT parameters list ex.: [['Na_+1', 'Cl_-1', 0.03]]
    """
    substances = None
    charges = None
    cations = None
    anions = None
    epsilon_matrix = None

    def __init__(self, ph: dict, parameters=sit_parameters):
        self.substances = np.array(ph)
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

    def get_y(self, ph: dict[str, float], T=298) -> dict[str, float]:
            """Calculate activity coefficients

            Parameters
            ----------
            ph: dict[str, float]
                Input dictionary
            T: float
                Temperature, [K], defaults to 298
            Returns
            -------
            dict[str, float]
                activity coefficients
            """
            molalities = np.array(list(ph.values()))
            I = get_I(molalities, charges=self.charges)
            sqrt_I = math.sqrt(I)
            D = 0.509 * sqrt_I / (1 + 1.5 * sqrt_I)

            lny = np.zeros((len(self.substances)))
            for i in range(len(self.substances)):
                if self.charges[i] > 0:
                    j = np.where(self.cations == self.substances[i])[0][0]
                    lny[i] = - self.charges[i] ** 2 * D + self.epsilon_matrix[j, :] @ molalities[self.charges < 0]
                elif self.charges[i] < 0:
                    j = np.where(self.anions == self.substances[i])[0][0]
                    lny[i] = - self.charges[i] ** 2 * D + self.epsilon_matrix[:, j] @ molalities[self.charges > 0]
                else:
                    lny[i] = 0
            y = {}
            for i in range(len(self.substances)):
                y[self.substances[i]] = 10 ** lny[i]
            return y

    def get_a(self, ph: dict[str, float], T=298):
        y = self.get_y(ph, T)
        a = []
        for i in range(len(self.substances)):
            a.append(ph[self.substances[i]] * y[self.substances[i]])
        return np.array(a)
