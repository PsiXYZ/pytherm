"""

This module contains a class :obj:`UNIFAC` for performing activity coefficient
calculations with the UNIFAC model.

.. contents:: :local:

UNIFAC Class
------------
.. autoclass:: UNIFAC
    :members: get_y, get_comb_original, get_comb_mod, get_res
    :undoc-members:
    :show-inheritance:
    :member-order: bysource


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
from math import log, exp, e
from typing import Callable

import numpy as np

from pytherm import constants
from .db import unifac as datasets
from .db.unifac import ParametersUNIFAC, SubstancesUNIFAC
from pytherm._activity import ActivityModel

R = constants.R


class UNIFAC(ActivityModel):
    r"""UNIFAC model for activity coefficients calculation

    UNIFAC type (classic or modified) depends on :obj:`.ParametersUNIFAC`

    Parameters
    ----------
    dataset : ParametersUNIFAC
        ParametersUNIFAC object with interaction parameters
    substances : SubstancesUNIFAC
        Substances UNIFAC object with substance's group representation

    Raises
    ------
    Warning
        Raises if unifac_mode is unknown (not classic or modified)

    Examples
    --------
    >>> import pytherm.activity.unifac as uf
    >>> subs = {
    ...    "n-hexane": "2*CH3 4*CH2",
    ...    "butanone-2": "1*CH3 1*CH2 1*CH3CO",
    ... }
    >>> system = {
    ...    'n-hexane': 0.5,
    ...    'butanone-2': 0.5,
    ... }
    >>> substances = uf.datasets.SubstancesUNIFAC()
    >>> substances.get_from_dict(subs)
    >>> am = uf.UNIFAC(dataset=uf.datasets.DOR, substances=substances)
    >>> am.get_y(conc=system, T=298)
    {'n-hexane': 1.514766775270851, 'butanone-2': 1.4331647782163541}
    """
    get_comb: Callable
    dict_mode: bool
    phase: SubstancesUNIFAC
    groups: list
    T = -1.0

    def __init__(self,
                 dataset: ParametersUNIFAC,
                 substances: SubstancesUNIFAC,
                 dict_mode=False):
        ActivityModel.__init__(self)
        
        self.unifac_mode = dataset['type']
        self.interaction_matrix = dataset['res']
        self.t_groups = dataset['comb']
        self.phase = substances
        self.dict_mode = dict_mode

        if self.unifac_mode == 'modified':
            self.get_comb = self.get_comb_mod
        elif self.unifac_mode == 'classic':
            self.get_comb = self.get_comb_original
        else:
            raise Warning("unknow unifac mode")

        self.__check_inter()
        self.__calculate_vdw()

        gr = []
        for i in substances:
            for j in substances[i].groups:
                if j not in gr:
                    gr.append(j)
        self.groups = gr

        gr_id = []
        for gr in self.groups:
            if self.t_groups[gr].id not in gr_id:
                gr_id.append(self.t_groups[gr].id)
        psi = {}
        for id1 in gr_id:
            b = {}
            if id1 not in psi:
                psi[id1] = {}
            for id2 in gr_id:
                b[id2] = 0
            psi[id1] = b
            self.psi = psi

            gamma_gr_pure = {}
            for comp in self.phase:
                gamma_gr_pure[comp] = {}
            for comp in self.phase:
                for gr in self.groups:
                    gamma_gr_pure[comp][gr] = 0
            self.gamma_gr_pure = gamma_gr_pure

    def get_y(self, system, T=298):
        r"""Calculate activity coefficietns for system dict

        .. math::
            \gamma_i =  \exp\left(\ln \gamma_i^c + \ln \gamma_i^r \right)

        Parameters
        ----------
        system
            Input dictionary {"Substance name": concentration}
        T : int, optional
            Temperature, [K], by default 298

        Returns
        -------
        dict[str, float]
            Activity coefficients
        """
        keys = list(self.phase.keys())

        if self.dict_mode:
            for i in system:
                self.phase[i].x = system[i]
        else:
            for i in range(len(self.phase)):
                self.phase[keys[i]].x = system[i]

        comb = self.get_comb(self.phase)
        res = self.get_res(self.phase, T)
        y = {}
        for i in keys:
            lny = comb[i] + res[i]
            y[i] = e ** lny

        if self.dict_mode:
            return y
        else:
            return np.array(list(y.values()))

    def get_comb_original(self, inp: SubstancesUNIFAC) -> dict[str, float]:
        r"""Calculate combinatorial component :math:`\ln\gamma_i^c` using original UNIFAC equation

        .. math::
            \ln \gamma_i^c = \ln \frac{\phi_i}{x_i} + \frac{z}{2} q_i
            \ln\frac{\theta_i}{\phi_i} + L_i - \frac{\phi_i}{x_i}
            \sum_{j=1}^{n} x_j L_j

        .. math::
            \theta_i = \frac{x_i q_i}{\sum_{j=1}^{n} x_j q_j}

        .. math::
            \phi_i = \frac{x_i r_i}{\sum_{j=1}^{n} x_j r_j}

        .. math::
            L_i = 5(r_i - q_i)-(r_i-1)

        Parameters
        ----------
        inp : SubstancesUNIFAC
            SubstancesUNIFAC object with defined concentrations
        z : int, optional
            Coordination number, by default 10

        Returns
        -------
        dict[str, float]
            Returns lny_c {"Substance name": lny_c}
        """

        phi = {}
        s = 0
        for comp in inp:
            s += self.phase[comp].r * self.phase[comp].x
        for comp in inp:
            phi[comp] = self.phase[comp].r / s

        theta = {}
        s = 0
        for comp in inp:
            s += self.phase[comp].q * self.phase[comp].x
        for comp in inp:
            theta[comp] = self.phase[comp].q / s

        rez = {}
        for i in inp:
            rez[i] = 1 - phi[i] + log(phi[i]) - 5 * self.phase[i].q * (
                        1 - phi[i] / theta[i] + log(phi[i] / theta[i]))
        return rez

    def get_comb_mod(self, inp: SubstancesUNIFAC) -> dict[str, float]:
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

        Parameters
        ----------
        inp : SubstancesUNIFAC
            SubstancesUNIFAC object with defined concentrations
        z : int, optional
            Coordination number, by default 10

        Returns
        -------
        dict[str, float]
            Returns lny_c {"Substance name": lny_c}
        """

        phi = {}
        s = 0
        for comp in inp:
            s += self.phase[comp].r * self.phase[comp].x
        for comp in inp:
            phi[comp] = self.phase[comp].r / s

        phi_m = {}
        s = 0
        for comp in inp:
            s += self.phase[comp].r ** (3/4) * self.phase[comp].x
        for comp in inp:
            phi_m[comp] = self.phase[comp].r ** (3/4) / s

        theta = {}
        s = 0
        for comp in inp:
            s += self.phase[comp].q * self.phase[comp].x
        for comp in inp:
            theta[comp] = self.phase[comp].q / s

        rez = {}
        for i in inp:
            rez[i] = 1 - phi_m[i] + log(phi_m[i]) - 5 * self.phase[i].q * (1 - phi[i]/theta[i] + log(phi[i]/theta[i]))
        return rez

    def get_res(self, inp, T: float):
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

        Parameters
        ----------
        inp : SubstancesUNIFAC
            SubstancesUNIFAC object with defined concentrations
        T : float
            Temperature, [K]

        Returns
        -------
        dict[str, float]
            Returns residual component {"Substance name": lny_a}
        """
        if self.T != T:
            self.calculate_psi(T)
            for comp in inp:
                x = inp[comp].x
                pure_comp = {
                    comp: inp[comp]
                }
                pure_comp[comp].x = 1
                self.gamma_gr_pure[comp] = self.__get_gamma_gr(pure_comp)
                pure_comp[comp].x = x
                self.T = T

        rez = {}
        gamma_gr = self.__get_gamma_gr(inp)
        for comp in inp:
            s = 0
            for gr in inp[comp].groups:
                s += inp[comp].groups[gr] * (gamma_gr[gr] - self.gamma_gr_pure[comp][gr])
            rez[comp] = s

        return rez

    def __get_gamma_gr(self, a):
        x = {}  # {имя группы: X}
        theta = {}  # {имя группы: тета}
        rez = {}

        s = 0
        for comp in a:
            for gr in a[comp].groups:
                s += a[comp].groups[gr] * a[comp].x
        for gr in self.groups:
            num = 0
            for comp in a:
                for gr2 in a[comp].groups:
                    if gr2 == gr:
                        num += a[comp].groups[gr2] * a[comp].x
            x[gr] = num / s

        # расчет тет
        s = 0
        for gr in self.groups:
            s += self.t_groups[gr].Q * x[gr]
        for gr in self.groups:
            theta[gr] = self.t_groups[gr].Q * x[gr] / s

        for gr_k in self.groups:
            s1 = 0
            # сумма1
            for gr_m in self.groups:
                m = self.t_groups[gr_m].id
                k = self.t_groups[gr_k].id
                psi_mk = self.psi[m][k]
                s1 += theta[gr_m] * psi_mk

            s2 = 0
            # сумма2
            for gr_m in self.groups:
                b = 0
                for gr_n in self.groups:
                    n = self.t_groups[gr_n].id
                    m = self.t_groups[gr_m].id
                    psi_nm = self.psi[n][m]
                    b += theta[gr_n] * psi_nm

                k = self.t_groups[gr_k].id
                m = self.t_groups[gr_m].id
                psi_km = self.psi[k][m]

                s2 += theta[gr_m] * psi_km / b

            rez[gr_k] = self.t_groups[gr_k].Q * (1 - log(s1) - s2)

        return rez

    def calculate_psi(self, T):
        for id1 in self.psi:
            for id2 in self.psi:
                a_ij = self.interaction_matrix[id1][id2]
                self.psi[id1][id2] = exp(- (a_ij[0]
                                          + a_ij[1] * T
                                          + a_ij[2] * T ** 2)
                                       / T)

    def get_ge(self, system: dict[str, float], T=298) -> float:
        """Calculate excess molar Gibbs free energy

        Parameters
        ----------
        system : dict[str, float]
            Input dictionary {"Substance name": concentration}
        T : int, optional
            Temperature, [K], by default 298

        Returns
        -------
        float
            excess molar Gibbs free energy
        """
        if self.dict_mode:
            y = self.get_y(system, T=T)
            ge = 0
            for sub in system:
                ge += system[sub] * log(y[sub])
            return R * T * ge
        else:
            y = self.get_y(system, T=T)
            ge = np.sum(system * np.log(y))
            return R * T * ge

    def get_ge_RT(self, system, T=298):
        """Calculate excess molar Gibbs free energy divided by RT

        Parameters
        ----------
        system : dict[str, float]
            Input dictionary {"Substance name": concentration}
        T : int, optional
            Temperature, [K], by default 298

        Returns
        -------
        float
            excess molar Gibbs free energy
        """
        if self.dict_mode:
            y = self.get_y(system, T=T)
            ge = 0
            for sub in system:
                ge += system[sub] * log(y[sub])
            return ge
        else:
            y = self.get_y(system, T=T)
            ge = np.sum(system * np.log(y))
            return ge

    def get_t2(self, i, j):
        print(self.interaction_matrix[i][j])

    def get_t1(self, name):
        print(self.t_groups[name].id,
              self.t_groups[name].R, self.t_groups[name].Q)

    def get_gr(self, name):
        return self.phase[name].groups

    def get_gr_str(self, name):
        s = ""
        g = self.phase[name].groups
        for i in g:
            s += f"{g[i]}*{i} "
        return (s[:-1])

    def __check_inter(self):
        """Checks for the presence of aij and aji parameters for all groups
        """
        groups = []
        ph = self.phase
        # создается список с именами всех групп в данной фазе
        for i in ph:
            for j in ph[i].groups:
                if j not in groups:
                    groups.append(j)

        for gr1 in groups:
            if self.t_groups[gr1].id not in self.interaction_matrix:
                self.__parameters_alert(gr1)
            for gr2 in groups:
                if self.t_groups[gr2].id not in self.interaction_matrix[self.t_groups[gr1].id]:
                    self.__parameters_alert(gr1, gr2)

    def __parameters_alert(self, *p):
        raise Warning(
            f"NO INTERACTION PARAMETERS FOR: {p}"
        )

    def change_t2(self, i, j, vals):
        self.interaction_matrix[i - 1][j - 1] = vals

    def __calculate_vdw(self):
        for comp in self.phase:
            r, q = 0, 0
            for gr in self.phase[comp].groups:
                r += self.phase[comp].groups[gr] * self.t_groups[gr].R
                q += self.phase[comp].groups[gr] * self.t_groups[gr].Q
            self.phase[comp].r = r
            self.phase[comp].q = q


class UNIFAC_W(UNIFAC):
    def __init__(
            self,
            dataset: ParametersUNIFAC,
            substances: SubstancesUNIFAC,
            Mw,
            dict_mode=False
    ):
        self.Mw = Mw
        super().__init__(dataset=dataset, substances=substances, dict_mode=dict_mode)

    def get_y(self, system, T=298.0):
        keys = list(system.keys())

        w_M = 0
        for i in keys:
            w_M += system[i] / self.Mw[i]

        system_x = {}
        for i in keys:
            system_x[i] = (system[i] / self.Mw[i]) / w_M

        y = super().get_y(system_x, T)

        y_w = {}
        for i in keys:
            y_w[i] = y[i] / (self.Mw[i] * w_M)

        return y_w
