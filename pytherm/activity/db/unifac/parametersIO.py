from pytherm.activity.db.unifac import defsubs


class ParametersUNIFAC(dict):
    def __init__(self, name):
        self['name'] = name

    def set_res(self, params: list, key="res"):
        if key not in self:
            self[key] = {}
        for i, j, *p in params:
            if i not in self[key]:
                self[key][i] = {}
            self[key][i][j] = p

    def set_comb(self, params: list, key="comb"):
        if key not in self:
            self[key] = {}
        for id, name, R, Q in params:
            self[key][name] = Group(id=id,
                                    r=R,
                                    q=Q)

    def set_type(self, type, key="type"):
        self[key] = type


class Group:
    """Class for unifac single group
    """

    def __init__(self, id: int, r: float, q: float):
        """

        Args:
            id (int): group id in parameter table
            r (float): volume parameter
            q (float): surface parameter
        """
        self.id = int(id)
        self.R = float(r)
        self.Q = float(q)


class Substance:
    """Class for unifac substance using Group
    """

    def __init__(self, concentration, groups):
        if concentration is not None:
            self.x = float(concentration)
        self.groups = groups  # словарь {Имя группы : кол-во}


class SubstancesUNIFAC(dict):
    from . import defsubs

    def __init__(self):
        pass

    def create_phase(self,
                     subs: dict,
                     ) -> dict[str, Substance]:
        """Create phase dict for Unifac

        Args:
            inp (dict): {substance_name: concentration}  dictionary

        Returns:
            dict[str, Substance]: dict for Unifac calculation
        """
        # словарь веществ - групп
        tsubs = self.get_subs(subs)
        for s in subs:
            gr = {}
            for i in tsubs[s]:
                gr[i[0]] = i[1]
            self[s] = Substance(0, gr)

    def get_subs(self, subs) -> dict:
        """

        Args:

        Returns:
            dict: {substance name: [Group, count]}
        """
        data = {}

        buf = subs
        for i in buf:
            b = buf[i].split(" ")
            a = []
            for j in b:
                # int -> float
                a.append([j.split("*")[1], float(j.split("*")[0])])
            data[i] = a
        return data

    def get_from_defsubs(self, subs):
        b = {}
        for s in subs:
            b[s] = defsubs.subs[s]
        self.create_phase(b)

    def get_from_dict(self, subs):
        b = {}
        for s in subs:
            b[s] = subs[s]
        self.create_phase(b)
           
    def get_from_csv(self,
                     subs,
                     path=r"D:\код\pytherm use\pytherm\pytherm\activity\db\unifac\sub.csv",
                    ):
        df = open(path, 'r')
        a = {}
        for line in df:
            buf = line[:-1].split(': ')
            a[buf[0]] = buf[1]
        b = {}
        for s in subs:
            b[s] = a[s]
        self.create_phase(b)
