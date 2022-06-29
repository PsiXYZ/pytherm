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
