from math import log, exp, e
from .db.unifac import io

R = 8.314462


class Unifac:
    def __init__(self, inp, mode="g"):
        self.mode = mode
        self.phase = self.create_ph(inp)  # словарь название в-ва - Sub из io

        self.inter = io.get_inter(mode)  # параметры взаимодействия
        self.t_groups = io.get_groups(mode)  # собственные параметры групп

    def get_y(self, inp, t=298):
        for i in inp:
            self.phase[i].x = inp[i]

        comb = self.get_comb(self.phase)
        res = self.get_res(self.phase, t)
        y = {}
        for i in inp:
            lny = comb[i] + res[i]
            y[i] = e ** lny

        return y

    def get_comb(self, inp, z=10):
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

    def get_res(self, inp, temp):
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
                    a_mk = self.inter[self.t_groups[t].id - 1][self.t_groups[s].id - 1]
                    s1 += tet[t] * exp(- a_mk / temp)

                # сумма под вторым логарифмом
                for t in gr:
                    a_km = self.inter[self.t_groups[s].id - 1][self.t_groups[t].id - 1]
                    num = tet[t] * exp(- a_km / temp)  # значение в числителе для m
                    den = 0
                    # расчет знаменателя
                    for u in gr:
                        a_nm = self.inter[self.t_groups[u].id - 1][self.t_groups[t].id - 1]
                        den += tet[u] * exp(- a_nm / temp)
                    s2 += num / den
                rez[s] = self.t_groups[s].Q * (1 - log(s1) - s2)
            return rez

        g = get_g(inp)

        # расчет Г в растворе содержащем только один тип молекул Гik
        for i in inp:
            buf = {}
            buf[i] = inp[i]  # создание фазы из одного компонента, концентрации можно не менять тк все нормируется
            gp[i] = get_g(buf)

        for i in inp:
            s = 0
            for k in inp[i].groups:
                s += inp[i].groups[k] * (g[k] - gp[i][k])
            rez[i] = s

        return rez

    def create_ph(self, inp):

        ph = io.create_ph(inp, self.mode)
        return ph

    def get_ge(self, inp, t=298, n=1):
        y = self.get_act(inp, t=t)
        ge = 0
        for sub in inp:
            ge += n*inp[sub]*log(y[sub])
        return R*t*ge

    def get_t2(self, i, j):
        print(self.inter[i - 1][j - 1])

    def get_t1(self, name):
        print(self.t_groups[name].id, self.t_groups[name].R, self.t_groups[name].Q)

    def get_gr(self, name):
        print(self.phase[name].groups)