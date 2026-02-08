"""

This module contains a class :obj:`UNIFAC` for performing activity coefficient
calculations with the UNIFAC model.

.. contents:: :local:

UNIFAC Class
------------
.. autoclass:: UNIFAC
    :members: get_y, get_y_array
    :undoc-members:
    :member-order: bysource

UNIFAC_W Class
--------------
.. autoclass:: UNIFAC_W
    :members: get_y
    :undoc-members:
    :member-order: bysource

Functions
---------
.. autofunction:: get_res
.. autofunction:: get_comb_classic
.. autofunction:: get_comb_mod

UNIFAC substances
-----------------
Substances must be a special :obj:`.SubstancesUNIFAC` object

UNIFAC parameters
-----------------
Parameters must be a special :obj:`.ParametersUNIFAC` object

There are some ready to use :obj:`.ParametersUNIFAC` objects in :obj:`.unifac.datasets`:

* Classic UNIFAC:
    * :obj:`.unifac.datasets.VLE` (:obj:`.db.unifac.parameters.vle.VLE`) [1]_
    * :obj:`.unifac.datasets.LLE` (:obj:`.db.unifac.parameters.vle.LLE`) [2]_
    * :obj:`.unifac.datasets.INF` (:obj:`.db.unifac.parameters.vle.INF`) [3]_
    * :obj:`.unifac.datasets.BIO2016_1` (:obj:`.db.unifac.parameters.vle.BIO2016_1`) [6]_
    * :obj:`.unifac.datasets.BIO2016_2` (:obj:`.db.unifac.parameters.vle.BIO2016_2`) [6]_
* Modified UNIFAC:
    * :obj:`.unifac.datasets.DOR` (:obj:`.db.unifac.parameters.vle.DOR`) [4]_
    * :obj:`.unifac.datasets.NIST2015` (:obj:`.db.unifac.parameters.vle.NIST2015`) [5]_


References
----------

.. [1] Published DDB parameters, 2021 JAN,
    https://www.ddbst.com/published-parameters-unifac.html
.. [2] Magnussen1981, DOI: https://doi.org/10.1021/i200013a024
.. [3] Bastos1988, DOI: https://doi.org/10.1021/i200013a024
.. [4] DOR, published DDB parameters, 2021 JAN,
    https://www.ddbst.com/PublishedParametersUNIFACDO.html
.. [5] Kang2015, DOI: https://doi.org/10.1016/j.fluid.2014.12.042
.. [6] Bessa2016, DOI: https://doi.org/10.1016/j.fluid.2016.05.020

"""
from numba import njit
import numpy as np
from pytherm import constants
from .db import unifac as datasets
from .db.unifac import ParametersUNIFAC, SubstancesUNIFAC
from pytherm.cpp import ActivityModel

R = constants.R
use_numba_cache = True


class UNIFAC(ActivityModel):
    r"""UNIFAC model for activity coefficients calculation

        UNIFAC type (classic or modified) depends on :obj:`.ParametersUNIFAC`
        For classic mode UNIFAC use :func:`get_comb_classic`, for modified :func:`get_comb_mod`

        Parameters
        ----------
        dataset : ParametersUNIFAC
            ParametersUNIFAC object with interaction parameters
        substances : SubstancesUNIFAC
            Substances UNIFAC object with substance's group representation

        Examples
        --------
        >>> import pytherm.activity.unifac as uf
        >>> subs = {
        ...    "n-hexane": "2*CH3 4*CH2",
        ...    "butanone-2": "1*CH3 1*CH2 1*CH3CO",
        ... }
        >>> x = [0.5, 0.5]
        >>> substances = uf.datasets.SubstancesUNIFAC()
        >>> substances.get_from_dict(subs)
        >>> am = uf.UNIFAC(dataset=uf.datasets.DOR, substances=substances)
        >>> am.get_y(x=x, T=298)
        {'n-hexane': 1.514766775270851, 'butanone-2': 1.4331647782163541}
        """

    gr_names: np.ndarray  # array with group names
    gr_id_global: np.ndarray  # array with global main id for each group from ParametersUNIFAC
    gr_id_local: np.ndarray  # array with local main id for each group (from 0 to 1)
    group_comp: np.ndarray[(np.any, np.any,)]  # matrix with group composition
    matrix_id: np.ndarray  # array with global gr main id
    res_matrix: np.ndarray[(np.any, np.any,), np.any]  # matrix with a, b, c for res component
    group_R: np.ndarray  # R array (for group)
    group_Q: np.ndarray  # Q array (for group)
    comp_r: np.ndarray  # q array (for component)
    comp_q: np.ndarray  # r array (for component)
    n_gr: int  # number of UNIFAC groups
    n_comp: int  # number of components
    T: float  # current temperature, K
    psi: np.ndarray[(np.any, np.any,)]  # psi matrix for current temperature
    ln_gamma_pure: np.ndarray[(np.any, np.any,)]  # matrix of groups ln_gamma for pure components
    modified_mode: bool  # if True use modified combinatorial part

    def __init__(self, dataset: ParametersUNIFAC, substances: SubstancesUNIFAC):
        ActivityModel.__init__(self)
        
        gr_list = []
        id_list = []
        for s in substances:
            for gr in substances[s].groups:
                if gr not in gr_list:
                    gr_list.append(gr)
                    id_list.append(dataset['comb'][gr].id)
        self.gr_names = np.array(gr_list)
        self.gr_id_global = np.array(id_list)

        comp_list = []
        for s in substances:
            b = []
            for gr in self.gr_names:
                if gr in substances[s].groups:
                    b.append(substances[s].groups[gr])
                else:
                    b.append(0)
            comp_list.append(b)
        self.group_comp = np.array(comp_list)

        q_list = []
        r_list = []
        for gr in self.gr_names:
            q_list.append(dataset['comb'][gr].Q)
            r_list.append(dataset['comb'][gr].R)
        self.group_Q = np.array(q_list)
        self.group_R = np.array(r_list)

        id = []
        for i in self.gr_id_global:
            if i not in id:
                id.append(i)
        self.matrix_id = np.array(id)

        loc_id = []
        for i, gr_i in enumerate(self.gr_id_global):
            for j, gr_j in enumerate(self.matrix_id):
                if gr_i == gr_j:
                    loc_id.append(j)
        self.gr_id_local = np.array(loc_id)

        inter = []
        for i in self.matrix_id:
            b = []
            for j in self.matrix_id:
                if j in dataset['res'][i]:
                    b.append(dataset['res'][i][j])
                else:
                    self.__parameters_alert(i, j)
            inter.append(b)
        self.res_matrix = np.array(inter)

        self.n_gr = len(self.gr_id_global)
        self.n_comp = len(substances)

        self.comp_r, self.comp_q = self.__get_vdw_params()

        self.T = -1

        n_main_gr = len(self.matrix_id)
        self.psi = np.zeros((n_main_gr, n_main_gr))
        self.ln_gamma_pure = np.zeros((self.n_comp, self.n_gr))

        self.unifac_mode = dataset['type']
        if self.unifac_mode == 'modified':
            self.modified_mode = True
        elif self.unifac_mode == 'classic':
            self.modified_mode = False
        else:
            raise Warning("unknown unifac mode")

    def get_y(self, conc: np.ndarray, T=298.0) -> np.ndarray:
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
        >>> UNIFAC.get_y([0.5, 0.5], T=298)
        """
        if self.T != T:
            self.T = T
            self.psi = get_psi(T, self.matrix_id, self.res_matrix)
            self.ln_gamma_pure = get_gamma_pure(conc, self.n_comp, self.n_gr, self.group_comp, self.group_Q,
                                                self.gr_id_local, self.psi)
        y = get_y(conc, self.comp_r, self.comp_q, self.n_gr, self.group_comp, self.group_Q, self.gr_id_local,
                  self.n_comp, self.psi, self.ln_gamma_pure, self.modified_mode)
        return y

    def get_y_array(self, conc: np.ndarray[(np.any, np.any,)], T=298) -> np.ndarray[(np.any, np.any,)]:
        r"""Calculate activity coefficients for conc matrix

        .. math::
            \gamma_i =  \exp\left(\ln \gamma_i^c + \ln \gamma_i^r \right)

        Parameters
        ----------
        conc : np.ndarray
            Input concentration matrix, [molar fraction]
        T : float, optional
            Temperature, [K], by default 298.0

        Returns
        -------
        np.ndarray[(np.any, np.any,)]
            Activity coefficients

        Examples
        --------
        >>> UNIFAC.get_y([[0.5, 0.5], [[0.6, 0.4]]])
        """
        if self.T != T:
            self.T = T
            self.psi = get_psi(T, self.matrix_id, self.res_matrix)
            self.ln_gamma_pure = get_gamma_pure(conc, self.n_comp, self.n_gr, self.group_comp, self.group_Q,
                                                self.gr_id_local, self.psi)
        y = get_y_array(conc, self.comp_r, self.comp_q, self.n_gr, self.group_comp, self.group_Q, self.gr_id_local,
                        self.n_comp, self.psi, self.ln_gamma_pure, self.modified_mode)
        return y

    def get_ge(self, conc: np.ndarray, T=298.0) -> float:
        """Calculate excess molar Gibbs free energy

        Parameters
        ----------
        conc : np.ndarray
            Input concentration array, [molar fraction]
        T : int, optional
            Temperature, [K], by default 298

        Returns
        -------
        float
            Excess molar Gibbs free energy
        """
        y = self.get_y(conc, T)
        ge = np.sum(conc * np.log(y))
        return R * T * ge

    def __get_vdw_params(self) -> (np.ndarray, np.ndarray):
        """Calculate r and q parameters

        .. math::
            r_i = \sum_{k=1}^{n}\nu_kR_k
            q_i = \sum_{k=1}^{n}\nu_kQ_k


        Returns
        -------
        (np.ndarray, np.ndarray)
           (r array, q array)
        """
        r = np.zeros(self.n_comp)
        q = np.zeros(self.n_comp)
        for i in range(self.n_comp):
            r[i] = np.sum(self.group_R * self.group_comp[i])
            q[i] = np.sum(self.group_Q * self.group_comp[i])
        return r, q

    def __parameters_alert(self, *p):
        raise Warning(
            f"NO INTERACTION PARAMETERS FOR: {p}"
        )


class UNIFAC_W(UNIFAC):
    r"""UNIFAC model for activity coefficients calculation in weight fractions.

    Override :obj:`UNIFAC` methods to calculate properties in weight fractions

    Parameters
    ----------
    dataset : ParametersUNIFAC
        ParametersUNIFAC object with interaction parameters
    substances : SubstancesUNIFAC
        Substances UNIFAC object with substance's group representation
    molar_weight
        array of molar weighs

    """
    molar_weight: np.ndarray

    def __init__(
            self,
            dataset: ParametersUNIFAC,
            substances: SubstancesUNIFAC,
            molar_weight,
    ):
        self.molar_weight = molar_weight
        super().__init__(dataset=dataset, substances=substances)

    def get_y(self, conc: np.ndarray, T=298.0) -> np.ndarray:
        r"""Calculate activity coefficients for conc array

        Concentrations must be in weight fractions

        .. math::
            \gamma_i =  \exp\left(\ln \gamma_i^c + \ln \gamma_i^r \right)

        Parameters
        ----------
        conc : np.ndarray
            Input concentration array, [weight fraction]
        T : float, optional
            Temperature, [K], by default 298.0

        Returns
        -------
        np.ndarray
            Activity coefficients

        Examples
        --------
        """
        w_M = np.sum(conc / self.molar_weight)
        x = (conc / self.molar_weight) / w_M
        y = super().get_y(x, T)
        y_w = y / (self.molar_weight * w_M)
        return y_w


@njit(cache=use_numba_cache)
def get_y(conc: np.ndarray, comp_r: np.ndarray, comp_q: np.ndarray, n_gr: int,
          group_comp: np.ndarray[(np.any, np.any,)], group_Q: np.ndarray, gr_id_local: np.ndarray, n_comp: int,
          psi: np.ndarray[(np.any, np.any,)], ln_gamma_pure: np.ndarray[(np.any, np.any,)],
          modified_mode: bool) -> np.ndarray:
    if modified_mode:
        comb = get_comb_mod(conc, comp_r, comp_q)
    else:
        comb = get_comb_classic(conc, comp_r, comp_q)
    res = get_res(conc, n_gr, group_comp, group_Q, gr_id_local, n_comp, psi, ln_gamma_pure)

    lny = comb + res
    y = np.exp(lny)

    return y


@njit(cache=use_numba_cache)
def get_y_array(conc_array: np.ndarray[(np.any, np.any,)], comp_r: np.ndarray, comp_q: np.ndarray, n_gr: int,
                group_comp: np.ndarray[(np.any, np.any,)], group_Q: np.ndarray, gr_id_local: np.ndarray, n_comp: int,
                psi: np.ndarray[(np.any, np.any,)], ln_gamma_pure: np.ndarray[(np.any, np.any,)],
                modified_mode: bool):
    y = np.zeros((len(conc_array), 2))
    for i in range(len(conc_array)):
        if modified_mode:
            comb = get_comb_mod(conc_array[i], comp_r, comp_q)
        else:
            comb = get_comb_classic(conc_array[i], comp_r, comp_q)
        res = get_res(conc_array[i], n_gr, group_comp, group_Q, gr_id_local, n_comp, psi, ln_gamma_pure)
        lny = comb + res
        y[i] = np.exp(lny)

    return y


@njit(cache=use_numba_cache)
def get_comb_mod(conc: np.ndarray, comp_r: np.ndarray, comp_q: np.ndarray) -> np.ndarray:
    r"""Calculate combinatorial component :math:`\ln\gamma_i^c` using modified equation

    .. math::
        \ln \gamma_i^c = 1 - {V'}_i + \ln({V'}_i) - 5q_i \left(1
        - \frac{V_i}{F_i}+ \ln\left(\frac{V_i}{F_i}\right)\right)

    .. math::
        V'_i = \frac{r_i^{3/4}}{\sum_j r_j^{3/4}x_j}

    .. math::
        V_i = \frac{r_i}{\sum_j r_j x_j}

    .. math::
        F_i = \frac{q_i}{\sum_j q_j x_j}

    """
    phi_m = comp_r ** (3 / 4) / np.sum(comp_r ** (3 / 4) * conc)
    phi = comp_r / np.sum(comp_r * conc)
    theta = comp_q / np.sum(comp_q * conc)

    return 1 - phi_m + np.log(phi_m) - 5 * comp_q * (1 - phi / theta + np.log(phi / theta))


@njit(cache=use_numba_cache)
def get_comb_classic(conc: np.ndarray, comp_r: np.ndarray, comp_q: np.ndarray) -> np.ndarray:
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


@njit(cache=use_numba_cache)
def get_res(conc: np.ndarray, n_gr: int, group_comp: np.ndarray[(np.any, np.any,)], group_Q: np.ndarray,
            gr_id_local: np.ndarray, n_comp: int, psi: np.ndarray[(np.any, np.any,)],
            ln_gamma_pure: np.ndarray[(np.any, np.any,)]) -> np.ndarray:
    r"""Calculate residual component :math:`\ln\gamma_i^r` using original equation


    .. math::
        \ln \gamma_i^r = \sum_{k}^n \nu_k^{(i)} \left[ \ln \Gamma_k
        - \ln \Gamma_k^{(i)} \right]

    .. math::
        \ln \Gamma_k = Q_k \left[1 - \ln \sum_m \Theta_m \Psi_{mk} - \sum_m
        \frac{\Theta_m \Psi_{km}}{\sum_n \Theta_n \Psi_{nm}}\right]

    .. math::
        \Theta_m = \frac{Q_m X_m}{\sum_{n} Q_n X_n}

    .. math::
        X_m = \frac{ \sum_j \nu^j_m x_j}{\sum_j \sum_n \nu_n^j x_j}

    .. math::
        \Psi_{mn} = \exp\left(\frac{-a_{mn} - b_{mn}T - c_{mn}T^2}{T}\right)

    """
    ln_gamma = get_gamma_gr(conc, n_gr, group_comp, group_Q, gr_id_local, psi)
    y = np.zeros_like(conc)

    for i in range(n_comp):
        y[i] = np.sum(group_comp[i] * (ln_gamma - ln_gamma_pure[i]))
    return y


@njit(cache=use_numba_cache)
def get_gamma_gr(conc: np.ndarray, n_gr: int, group_comp: np.ndarray[(np.any, np.any,)], group_Q: np.ndarray,
                 gr_id_local: np.ndarray, psi: np.ndarray[(np.any, np.any,)]):
    s = np.array([np.sum(group_comp.T[i] * conc) for i in range(n_gr)])
    X = s / np.sum(s)
    theta = X * group_Q / np.sum(X * group_Q)

    # temps = np.array((1, T, T ** 2))
    # n_main_gr = len(matrix_id)
    # psi = np.zeros((n_main_gr, n_main_gr))
    # for i in range(n_main_gr):
    #     for j in range(n_main_gr):
    #         psi[i][j] = np.exp(np.sum(- res_matrix[i][j] * temps / T))

    gamma_gr = np.zeros_like(theta)
    for k in range(n_gr):
        s1 = 0
        for m in range(n_gr):
            s1 += theta[m] * psi[gr_id_local[m]][gr_id_local[k]]

        s2 = 0
        for m in range(n_gr):
            b = 0
            for n in range(n_gr):
                b += theta[n] * psi[gr_id_local[n]][gr_id_local[m]]
            s2 += theta[m] * psi[gr_id_local[k]][gr_id_local[m]] / b

        gamma_gr[k] = group_Q[k] * (1 - np.log(s1) - s2)
    return gamma_gr


@njit(cache=use_numba_cache)
def get_psi(T: float, matrix_id: np.ndarray,
            res_matrix: np.ndarray[(np.any, np.any,), np.any]) -> np.ndarray[(np.any, np.any,)]:
    temps = np.array((1, T, T ** 2))
    n_main_gr = len(matrix_id)
    psi = np.zeros((n_main_gr, n_main_gr))
    for i in range(n_main_gr):
        for j in range(n_main_gr):
            psi[i][j] = np.exp(np.sum(- res_matrix[i][j] * temps / T))
    return psi


@njit(cache=use_numba_cache)
def get_gamma_pure(conc: np.ndarray, n_comp: int, n_gr: int, group_comp: np.ndarray[(np.any, np.any,)],
                   group_Q: np.ndarray, gr_id_local: np.ndarray, psi: np.ndarray[(np.any, np.any,)]):
    ln_gamma_pure = np.zeros((n_comp, n_gr))
    for i in range(n_comp):
        pure_conc = np.zeros_like(conc)
        pure_conc[i] = 1
        ln_gamma_pure[i] = get_gamma_gr(pure_conc, n_gr, group_comp, group_Q, gr_id_local, psi)
    return ln_gamma_pure
