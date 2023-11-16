"""

This module contains a class :obj:`UNIQUAC` for performing activity coefficient calculations with the UNIQUAC model.

.. contents:: :local:

UNIQUAC Class
------------
.. autoclass:: UNIQUAC
    :members: get_y
    :undoc-members:
    :member-order: bysource

Functions
---------
.. autofunction:: get_res
.. autofunction:: get_comb

"""

from pytherm.cpp import ActivityModel
import numpy as np
from pytherm import constants
R = constants.R


class UNIQUAC(ActivityModel):
    r"""UNIFAC model for activity coefficients calculation


    Parameters
    ----------
    comp_r: np.ndarray
        Array with r
    comp_q: np.ndarray
        Array with q
    res_params: np.ndarray[(np.any, np.any,), np.any])
        Matrix with coefficients for res component

    Examples
    --------
    >>> import numpy as np
    >>> from pytherm.activity import uniquac as uq
    >>> T = 25.0 + 273.15
    >>> xs = [0.7273, 0.0909, 0.1818]
    >>> rs = [.92, 2.1055, 3.1878]
    >>> qs = [1.4, 1.972, 2.4]
    >>> inter = [
    ...     [[0, 0], [0, 526.02], [0, 309.64]],
    ...     [[0, -318.06], [0, 0], [0, -91.532]],
    ...     [[0, 1325.1], [0, 302.57], [0, 0]],
    ... ]
    >>> am = uq.UNIQUAC(rs, qs, inter)
    >>> am.get_y(np.array(xs), T)
    [ 1.57039333  0.29482416 18.11432905]
    """
    comp_r: np.ndarray  # r array
    comp_q: np.ndarray  # q array
    res_matrix: np.ndarray[(np.any, np.any,), np.any]  # matrix with coefficients for res component
    T: float  # current temperature, K
    n_comp: int  # number of components
    t_matrix: np.ndarray[(np.any, np.any,)]  # t matrix for current temperature

    def __init__(self, comp_r: np.ndarray, comp_q: np.ndarray, res_params: np.ndarray[(np.any, np.any,), np.any]):
        ActivityModel.__init__(self)
        
        self.comp_r = np.array(comp_r)
        self.comp_q = np.array(comp_q)
        self.T = -1
        self.n_comp = len(comp_r)
        self.res_matrix = - np.array(res_params)

    def get_y(self, conc: np.ndarray, T=298.0):
        r"""Calculate activity coefficients for conc array

        Concentrations must be in molar fractions

        .. math::
            \gamma_i =  \exp\left(\ln \gamma_i^c + \ln \gamma_i^r \right)

        Parameters
        ----------
        conc : np.ndarray
            Input concentration array, [molar fraction]
        T : float, optional
            Temperature, [K], by default 298.0

        Returns
        -------
        np.ndarray
            Activity coefficients

        Examples
        --------
        >>> UNIQUAC.get_y([0.5, 0.5], T=298)
        """
        if self.T != T:
            self.T = T
            self.t_matrix = get_tmatrix(T, self.res_matrix)
        return get_y(conc, self.comp_r, self.comp_q, self.t_matrix)


def get_y(conc, comp_r, comp_q, t_matrix):
    comb = get_comb(conc, comp_r, comp_q)
    res = get_res(conc, comp_q, t_matrix)
    lny = comb + res
    y = np.exp(lny)
    return y


def get_comb(conc: np.ndarray, comp_r: np.ndarray, comp_q: np.ndarray) -> np.ndarray:
    r"""Calculate combinatorial component :math:`\ln\gamma_i^c` using classic equation

    .. math::
        \ln \gamma_i^c = 1 - {V}_i + \ln({V}_i) - 5q_i \left(1
        - \frac{V_i}{F_i}+ \ln\left(\frac{V_i}{F_i}\right)\right)

    .. math::
        V_i = \frac{r_i}{\sum_j r_j x_j}

    .. math::
        F_i = \frac{q_i}{\sum_j q_j x_j}

    """
    phi = comp_r / np.sum(comp_r * conc)
    theta = comp_q / np.sum(comp_q * conc)

    return 1 - phi + np.log(phi) - 5 * comp_q * (1 - phi / theta + np.log(phi / theta))


def get_res(conc: np.ndarray, comp_q: np.ndarray, t_matrix: np.ndarray[(np.any, np.any,)]) -> np.ndarray:
    r"""Calculate residual component :math:`\ln\gamma_i^r`

    .. math::
        \tau_{ji} = \exp\left(\frac{-\Delta u_{ij}}{RT}\right)
    .. math::
        \ln \tau_{ij} =a_{ij}+\frac{b_{ij}}{T}+c_{ij}\ln T + d_{ij}T
        + \frac{e_{ij}}{T^2}
    .. math::
        \ln \gamma_i^{res} = q_i \left(1 - \ln\frac{\sum_j^N q_j x_j \tau_{ji}}
        {\sum_j^N q_j x_j}- \sum_j \frac{q_k x_j \tau_{ij}}{\sum_k q_k x_k
        \tau_{kj}}\right)


    """
    n = len(conc)
    ln_y_res = np.zeros(n)
    for i in range(n):
        s1 = np.sum(comp_q * conc * t_matrix[:, i]) / np.sum(comp_q * conc)

        s2 = .0
        for j in range(n):
            s2 += comp_q[j] * conc[j] * t_matrix[i, j] / np.sum(comp_q * conc * t_matrix[:, j])
        ln_y_res[i] = comp_q[i] * (1 - np.log(s1) - s2)
    return ln_y_res


def get_tmatrix(T, res_matrix):
    temp = np.array([1, 1/T])
    n_components = len(res_matrix)
    t_matrix = np.zeros((n_components, n_components))
    for i in range(n_components):
        for j in range(n_components):
            s = np.sum(res_matrix[i][j] * temp)
            t_matrix[i][j] = np.exp(s)
    return t_matrix
