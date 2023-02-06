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
from .activity_model import ActivityModel

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
    >>> am.get_y(system=system, T=298)
    {'n-hexane': 1.514766775270851, 'butanone-2': 1.4331647782163541}
    """
    get_comb: Callable
    dict_mode: bool

    def __init__(self,
                 dataset: ParametersUNIFAC,
                 substances: SubstancesUNIFAC,
                 dict_mode=False):
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

    def get_y(self,
              system,
              T=298) -> dict[str, float]:
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

    def get_comb_original(self, inp: SubstancesUNIFAC, z=10) -> dict[str, float]:
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

        def get_r(name):
            r = 0
            for i in inp[name].groups:
                r += self.t_groups[i].R * inp[name].groups[i]
            return r

        def get_q(name):
            q = 0
            for i in inp[name].groups:
                q += self.t_groups[i].Q * inp[name].groups[i]
            return q

        def get_l(name):
            return z / 2 * (get_r(name) - get_q(name)) - (get_r(name) - 1)

        def get_f(name):
            s = 0
            for j in inp:
                rj = get_r(j)
                s += rj * inp[j].x
            ri = get_r(name)
            return ri * inp[name].x / s

        def get_t(name):
            s = 0
            for j in inp:
                qj = get_q(j)
                s += qj * inp[j].x
            qi = get_q(name)
            return qi * inp[name].x / s

        rez = {}
        for i in inp:
            qi = get_q(i)
            li = get_l(i)
            fi = get_f(i)
            ti = get_t(i)
            xi = inp[i].x

            s = 0
            for j in inp:
                s += inp[j].x * get_l(j)

            lny = log(fi / xi) + z / 2 * qi * log(ti / fi) + li - fi / xi * s
            rez[i] = lny
        return rez

    def get_comb_mod(self, inp: SubstancesUNIFAC, z=10) -> dict[str, float]:
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

        def get_r(name):
            r = 0
            for i in inp[name].groups:
                r += self.t_groups[i].R * inp[name].groups[i]
            return r

        def get_q(name):
            q = 0
            for i in inp[name].groups:
                q += self.t_groups[i].Q * inp[name].groups[i]
            return q

        def get_fc(name):
            s = 0
            for j in inp:
                rj = get_r(j)
                s += (rj ** (3 / 4)) * inp[j].x
            ri = get_r(name)
            return (ri ** (3 / 4)) / s

        def get_f(name):
            s = 0
            for j in inp:
                rj = get_r(j)
                s += rj * inp[j].x
            ri = get_r(name)
            return ri * inp[name].x / s

        def get_t(name):
            s = 0
            for j in inp:
                qj = get_q(j)
                s += qj * inp[j].x
            qi = get_q(name)
            return qi * inp[name].x / s

        rez = {}
        for i in inp:
            qi = get_q(i)
            fi = get_f(i)
            fic = get_fc(i)
            ti = get_t(i)
            xi = inp[i].x

            lny = 1 - fic + log(fic) - z / 2 * qi * \
                  (log(fi / ti) + 1 - fi / ti)
            rez[i] = lny
        return rez

    def get_res(self, inp, T: float) -> dict[str, float]:
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
        g = {}  # группа: Г
        gp = {}  # вещество: {группа: Гi}
        rez = {}

        # phase -> {группа: Гi}
        def get_g(a):
            x = {}  # {имя группы: X}
            tet = {}  # {имя группы: тета}
            gr = []
            rez = {}

            # создается список с именами всех групп в данной фазе
            for i in a:
                for j in a[i].groups:
                    if j not in gr:
                        gr.append(j)

            # расчет х
            den = 0
            # расчет знаменателя
            for i in a:
                for j in a[i].groups:
                    den += a[i].groups[j] * a[i].x
            # расчет числителей
            for g in gr:
                num = 0
                for i in a:
                    for j in a[i].groups:
                        if j == g:
                            num += a[i].groups[j] * a[i].x
                x[g] = num / den

            # расчет тет
            den = 0
            # расчет знаменателя
            for t in gr:
                den += self.t_groups[t].Q * x[t]
            for t in gr:
                tet[t] = self.t_groups[t].Q * x[t] / den

            for s in gr:
                s1 = 0
                s2 = 0
                # сумма под первым логарифмом
                for t in gr:
                    m = self.t_groups[t].id
                    k = self.t_groups[s].id
                    a_mk = self.interaction_matrix[m][k]
                    s1 += tet[t] * exp(- (a_mk[0]
                                          + a_mk[1] * T
                                          + a_mk[2] * T ** 2)
                                       / T)

                # сумма под вторым логарифмом
                for t in gr:
                    k = self.t_groups[s].id
                    m = self.t_groups[t].id
                    a_km = self.interaction_matrix[k][m]
                    num = tet[t] * exp(
                        - (a_km[0]
                           + a_km[1] * T
                           + a_km[2] * T ** 2)
                        / T)  # значение в числителе для m
                    den = 0
                    # расчет знаменателя
                    for u in gr:
                        n = self.t_groups[u].id
                        m = self.t_groups[t].id
                        a_nm = self.interaction_matrix[n][m]
                        den += tet[u] * exp(- (a_nm[0]
                                               + a_nm[1] * T
                                               + a_nm[2] * T ** 2)
                                            / T)
                    s2 += num / den
                rez[s] = self.t_groups[s].Q * (1 - log(s1) - s2)
            return rez

        g = get_g(inp)

        # расчет Г в растворе содержащем только один тип молекул Гik
        for i in inp:
            buf = {}
            '''
            создание фазы из одного компонента,
            концентрации можно не менять тк все нормируется
            '''
            buf[i] = inp[i]
            gp[i] = get_g(buf)

        for i in inp:
            s = 0
            for k in inp[i].groups:
                s += inp[i].groups[k] * (g[k] - gp[i][k])
            rez[i] = s

        return rez

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
        y = self.get_y(system, T=T)
        ge = 0
        for sub in system:
            ge += system[sub] * log(y[sub])
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
        y = self.get_y(system, T=T)
        ge = 0
        for sub in system:
            ge += system[sub] * log(y[sub])
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
