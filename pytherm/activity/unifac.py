"""
parameter sources:
    VLE: published DDB parameters, 2021 JAN,
        https://www.ddbst.com/published-parameters-unifac.html
    LLE: Magnussen1981,
        DOI: https://doi.org/10.1021/i200013a024
    INF: Bastos1988,
        DOI: https://doi.org/10.1021/ie00079a030
    DOR: published DDB parameters, 2021 JAN,
        https://www.ddbst.com/PublishedParametersUNIFACDO.html
    NIST2015: Kang2015,
        DOI: https://doi.org/10.1016/j.fluid.2014.12.042

    VLE, LLE, INF: original COMB and RES (1 param)
    DOR, NIST2015: mod COMB and res (3 param)
"""
from math import log, exp, e

from pytherm import constants
from .db import unifac as datasets

R = constants.R


class Unifac:
    """Unifac model for activity coeffcient calculation
    """
    get_comb = None

    def __init__(self,
                 dataset: datasets.ParametersUNIFAC,
                 substances: datasets.SubstancesUNIFAC):

        self.unifac_mode = dataset['type']
        self.interaction_matrix = dataset['res']
        self.t_groups = dataset['comb']
        self.phase = substances

        if self.unifac_mode == 'modified':
            self.get_comb = self.get_comb_mod
        elif self.unifac_mode == 'classic':
            self.get_comb = self.get_comb_original
        else:
            raise Warning("unknow unifac mode")

        self.__check_inter()

    def get_y(self,
              inp: dict[str, float],
              temperature=298) -> dict[str, float]:
        """Calculate activity coefficient for input dict

        Args:
            inp (dict[str, float]): input dictionary {name: conentration}
            temperature (optional): Temperature in K. Defaults to 298.

        Returns:
            dict[str, float]: activity coefficient for inp substances
        """
        for i in inp:
            self.phase[i].x = inp[i]

        comb = self.get_comb(self.phase)
        res = self.get_res(self.phase, temperature)
        y = {}
        for i in inp:
            # print(i, e ** comb[i], e ** res[i])
            lny = comb[i] + res[i]
            y[i] = e ** lny
        return y

    def get_comb_original(self, inp: dict[str, float], z=10) -> dict:
        """Calculate combinatorial component using original equation

        Args:
            inp (dict[str, float]): input dictionary {name: conentration}
            z (optional): Coordination number. Defaults to 10.

        Returns:
            dict: combinatorial component
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

    def get_comb_mod(self, inp, z=10) -> dict:
        """Calculate combinatorial component using modified equation

        Args:
            inp (dict[str, float]): input dictionary {name: conentration}
            z (optional): Coordination number. Defaults to 10.

        Returns:
            dict: combinatorial component
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

    def get_res(self, inp, temperature) -> dict:
        """Calculate combinatorial component using original equation

        Args:
            inp (dict[str, float]): input dictionary {name: conentration}
            temperature (float): Temperature in K.

        Returns:
            dict: resudual component
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
                                          + a_mk[1] * temperature
                                          + a_mk[2] * temperature ** 2)
                                       / temperature)

                # сумма под вторым логарифмом
                for t in gr:
                    k = self.t_groups[s].id
                    m = self.t_groups[t].id
                    a_km = self.interaction_matrix[k][m]
                    num = tet[t] * exp(
                        - (a_km[0]
                           + a_km[1] * temperature
                           + a_km[2] * temperature ** 2)
                        / temperature)  # значение в числителе для m
                    den = 0
                    # расчет знаменателя
                    for u in gr:
                        n = self.t_groups[u].id
                        m = self.t_groups[t].id
                        a_nm = self.interaction_matrix[n][m]
                        den += tet[u] * exp(- (a_nm[0]
                                               + a_nm[1] * temperature
                                               + a_nm[2] * temperature ** 2)
                                            / temperature)
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

    def get_ge(self, inp, t=298, n=1):
        y = self.get_y(inp, temperature=t)
        ge = 0
        for sub in inp:
            ge += n*inp[sub]*log(y[sub])
        return R*t*ge

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
