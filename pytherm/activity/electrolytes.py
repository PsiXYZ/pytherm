import numpy as np


def get_I(molalities: np.ndarray, charges: np.ndarray):
    I = molalities @ charges ** 2
    return I / 2


def get_A(self, T=298):
    return (0.13422 *
            (
                4.1725332
                - 0.1481291 * T ** (0.5)
                + 1.5188505 * 10 ** (-5) * T ** 2
                - 1.8016317 * 10 ** (-8) * T ** 3
                + 9.3816144 * 10 ** (-10) * T ** (3.5)
            )
            )
