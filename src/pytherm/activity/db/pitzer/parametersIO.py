import numpy as np


class ParametersPitzer(dict):
    def __init__(self, parameters_name):
        self["parameters_name"] = parameters_name

    def set_ca(self, param: list, key="ca"):
        if key not in self:
            self[key] = {}
        for c, a, *ca in param:
            if c not in self[key]:
                self[key][c] = {}
            self[key][c][a] = ca

    def set_cc(self, param: list, key="cc"):
        if key not in self:
            self[key] = {}
        for c1, c2, *cc in param:
            if c1 not in self[key]:
                self[key][c1] = {}
            self[key][c1][c2] = cc[0]

    def set_aa(self, param: list, key="aa"):
        if key not in self:
            self[key] = {}
        for a1, a2, *aa in param:
            if a1 not in self[key]:
                self[key][a1] = {}
            self[key][a1][a2] = aa[0]

    def set_cca(self, param: list, key="cca"):
        if key not in self:
            self[key] = {}
        for c1, c2, a, *cca in param:
            if c1 not in self[key]:
                self[key][c1] = {}
            if c2 not in self[key][c1]:
                self[key][c1][c2] = {}
            self[key][c1][c2][a] = cca[0]

    def set_caa(self, param: list, key="caa"):
        if key not in self:
            self[key] = {}
        for c, a1, a2, *caa in param:
            if c not in self[key]:
                self[key][c] = {}
            if a1 not in self[key][c]:
                self[key][c][a1] = {}
            self[key][c][a1][a2] = caa[0]

    def set_nc(self, param: list, key="nc"):
        if key not in self:
            self[key] = {}
        for n, c, *nc in param:
            if n not in self[key]:
                self[key][n] = {}
            self[key][n][c] = nc[0]

    def set_na(self, param: list, key="na"):
        if key not in self:
            self[key] = {}
        for n, a, *na in param:
            if n not in self[key]:
                self[key][n] = {}
            self[key][n][a] = na[0]

    def set_nca(self, param: list, key="nca"):
        if key not in self:
            self[key] = {}
        for n, c, a, *nca in param:
            if n not in self[key]:
                self[key][n] = {}
            if c not in self[key][n]:
                self[key][n][c] = {}
            self[key][n][c][a] = nca[0]

    def set_nn(self, param: list, key="nn"):
        if key not in self:
            self[key] = {}
        for n1, n2, *nn in param:
            if n1 not in self[key]:
                self[key][n1] = {}
            self[key][n1][n2] = nn[0]

    def set_nnn(self, param: list, key="nnn"):
        if key not in self:
            self[key] = {}
        for n1, n2, n3, *nnn in param:
            if n1 not in self[key]:
                self[key][n1] = {}
            if n2 not in self[key][n1]:
                self[key][n1][n2] = {}
            self[key][n1][n2][n3] = nnn[0]

    def is_exist(self, key, *args):
        if len(args) == 1:
            if args[0] in self[key]:
                return True
            else:
                False

        if len(args) == 3:
            if args[0] in self[key]:
                if args[1] in self[key][args[0]]:
                    if args[2] in self[key][args[0]][args[1]]:
                        return True
                    else:
                        return False
                else:
                    return False
            else:
                return False

        if len(args) == 2:
            if args[0] in self[key]:
                if args[1] in self[key][args[0]]:
                    return True
                else:
                    return False
            else:
                return False

        if len(args) == 3:
            if args[0] in self[key]:
                if args[1] in self[key][args[0]]:
                    if args[2] in self[key][args[0]][args[1]]:
                        return True
                    else:
                        return False
                else:
                    return False
            else:
                return False


class ParametersPitzerNew(dict):
    k = {
        'B0': 2,
        'B1': 2,
        'B2': 2,
        'C0': 2,
        'THETA': 2,
        'LAMDA': 2,
        'ZETA': 3,
        'PSI': 3,
    }
    raw_parameters = {}
    poly_form = None

    def load_from_dat(self, dat_path, f):
        self.poly_form = f

        df = open(dat_path)
        bufer = []
        for line in df:
            b = line
            while '\t' in b:
                b = b.replace('\t', ' ')
            while '  ' in b:
                b = b.replace('  ', ' ')
            if b[0] == ' ':
                b = b[1:]
            if b[0] == '#':
                continue
            if '#' in b:
                i = b.index('#')
                b = b[:i]
            if b[-1:] == '\n':
                b = b[:-1]
            if b[-1] == ' ':
                b = b[:-1]
            bufer.append(b)
        df.close()

        for line in bufer:
            if line[0] == '-':
                key = line[1:]
                continue
            self.__unbox_values(line.split(' '), key)

        for key in self.raw_parameters:
            if key not in ['B0', 'B1', 'B2', 'C0']:
                n_subs = self.k[key]
                if n_subs == 2:
                    if key not in self:
                        self[key] = {}
                    self[key] = {s: {} for s in self.raw_parameters[key]}
                elif n_subs == 3:
                    if key not in self:
                        self[key] = {}
                    for s1 in self.raw_parameters[key]:
                        if s1 not in self[key]:
                            self[key][s1] = {}
                            for s2 in self.raw_parameters[key][s1]:
                                if s2 not in self[key][s1]:
                                    self[key][s1][s2] = {}
        self['CA'] = {i: {} for i in self.raw_parameters['B0']}
        for s1 in self.raw_parameters['B0']:
            for s2 in self.raw_parameters['B0'][s1]:
                self['CA'][s1][s2] = np.zeros((4), dtype=np.float64)

    def __unbox_values(self, line, key):
        n_subs = self.k[key]
        values = line[n_subs:]

        if n_subs == 2:
            s1, s2 = line[:n_subs]
            if key not in self.raw_parameters:
                self.raw_parameters[key] = {}
            if s1 not in self.raw_parameters[key]:
                self.raw_parameters[key][s1] = {}
            self.raw_parameters[key][s1][s2] = np.array(list(map(float, values)))
        elif n_subs == 3:
            s1, s2, s3 = line[:n_subs]
            if key not in self.raw_parameters:
                self.raw_parameters[key] = {}
            if s1 not in self.raw_parameters[key]:
                self.raw_parameters[key][s1] = {}
            if s2 not in self.raw_parameters[key][s1]:
                self.raw_parameters[key][s1][s2] = {}
            self.raw_parameters[key][s1][s2][s3] = np.array(list(map(float, values)))

    def set_poly(self, f):
        self.poly_form = f

    def update_params(self, T):
        coefs = self.poly_form(T)
        for key in self.raw_parameters:
            if key not in ['B0', 'B1', 'B2', 'C0']:
                n_subs = self.k[key]
                if n_subs == 2:
                    for s1 in self.raw_parameters[key]:
                        for s2 in self.raw_parameters[key][s1]:
                            params = self.raw_parameters[key][s1][s2]
                            poly = coefs[:len(params)]
                            self[key][s1][s2] = params @ poly
                elif n_subs == 3:
                    for s1 in self.raw_parameters[key]:
                        for s2 in self.raw_parameters[key][s1]:
                            for s3 in self.raw_parameters[key][s1][s2]:
                                params = self.raw_parameters[key][s1][s2][s3]
                                poly = coefs[:len(params)]
                                self[key][s1][s2][s3] = params @ poly

        p = ['B0', 'B1', 'B2', 'C0']
        for s1 in self.raw_parameters['B0']:
            for s2 in self.raw_parameters['B0'][s1]:
                for i in range(4):
                    key = p[i]
                    if s1 in self.raw_parameters[key]:
                        if s2 in self.raw_parameters[key][s1]:
                            params = self.raw_parameters[key][s1][s2]
                            poly = coefs[:len(params)]
                            self['CA'][s1][s2][i] = params @ poly

    def is_exist(self, key, *args):
        if len(args) == 1:
            if args[0] in self[key]:
                return True
            else:
                False

        if len(args) == 3:
            if args[0] in self[key]:
                if args[1] in self[key][args[0]]:
                    if args[2] in self[key][args[0]][args[1]]:
                        return True
                    else:
                        return False
                else:
                    return False
            else:
                return False

        if len(args) == 2:
            if args[0] in self[key]:
                if args[1] in self[key][args[0]]:
                    return True
                else:
                    return False
            else:
                return False

        if len(args) == 3:
            if args[0] in self[key]:
                if args[1] in self[key][args[0]]:
                    if args[2] in self[key][args[0]][args[1]]:
                        return True
                    else:
                        return False
                else:
                    return False
            else:
                return False
