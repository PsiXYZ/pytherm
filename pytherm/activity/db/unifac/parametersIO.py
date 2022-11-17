"""Subpackage for UNIFAC parameters and substance operation

.. contents:: :local:

Parameters processing
---------------------

.. autoclass:: ParametersUNIFAC
    :members: set_res, set_comb, set_type
    :show-inheritance:
    :member-order: bysource

.. autoclass:: Group
    :members:
    :undoc-members:
    :show-inheritance:
    :member-order: bysource

Substances processing
---------------------

.. autoclass:: Substance
    :members:
    :undoc-members:
    :show-inheritance:
    :member-order: bysource

.. autoclass:: SubstancesUNIFAC
    :members: get_from_dict
    :undoc-members:
    :show-inheritance:
    :member-order: bysource

"""

from pytherm.activity.db.unifac import defsubs


class ParametersUNIFAC(dict):
    """Class for UNIFAC parameters processing

    Res parameters are stored as dict['res'][i][j]. Comb parameters are stored as dict['comb']['name']: :obj:`Group`

    Examples
    --------
    >>> t1 = [
    ...    [1, 'CH3', 0.9011, 0.848],
    ...    [1, 'CH2', 0.6744, 0.54],
    ...    [1, 'CH', 0.4469, 0.228],
    ...    [1, 'C', 0.2195, 0.0],
    ...    [3, 'ACH', 0.5313, 0.4],
    ...    [3, 'AC', 0.3652, 0.12],
    ... ]
    >>> t2 = [
    ...    [1, 3, 61.13, 0.0, 0.0],
    ...    [3, 1, -11.12, 0.0, 0.0]
    >>> VLE = ParametersUNIFAC("VLE")
    >>> VLE.set_type("classic")
    >>> VLE.set_comb(t1)
    >>> VLE.set_res(t2)
    """

    def __init__(self, name):
        self['name'] = name

    def set_res(self, params: list, key="res"):
        r"""Set parameters for :math:`\ln\gamma_i^r` calculations

        If one of parameters is absent, then it must be zero.

        Parameters
        ----------
        params : list
            Parameters list [i, j, aij, bij, cij]
        """
        if key not in self:
            self[key] = {}
        for i, j, *p in params:
            if i not in self[key]:
                self[key][i] = {}
            self[key][i][j] = p

    def set_comb(self, params: list, key="comb"):
        r"""Set parameters for :math:`\ln\gamma_i^c` calculations

        Parameters
        ----------
        params : list
            Parameters list [id, 'group name', r, q],
        """
        if key not in self:
            self[key] = {}
        for id, name, R, Q in params:
            self[key][name] = Group(id=id,
                                    r=R,
                                    q=Q)

    def set_type(self, type, key="type"):
        """Set model type (classic or modified)

        Parameters
        ----------
        type : str
            classic or modified
        """
        self[key] = type


class Group:
    """Class for UNIFAC single group

    Contain group id, surface and volume

    Parameters
    ----------
    id : int
        group id in parameters table
    r : float
        group volume
    q : float
        group surface

    Attributes
    ----------
    id : int
        group id in parameters table
    r : float
        group volume
    q : float
        group surface
    """

    def __init__(self, id: int, r: float, q: float):
        self.id = int(id)
        self.R = float(r)
        self.Q = float(q)


class Substance:
    """Class for UNIFAC single substance

    Parameters
    ----------
    x : float
        concentration
    groups : dict[str, float]
        {'group name': count}

    Attributes
    ----------
    x : float
        Concentration
    groups : dict[str, float]
        {'group name': count}

    Examples
    --------
    >>> hexane = Substance(1, {'CH3': 2, 'CH2': 4})
    """

    def __init__(self, x: float, groups: dict[str, float]):

        if x is not None:
            self.x = float(x)
        else:
            self.x = 0
        self.groups = groups  # словарь {Имя группы : кол-во}


class SubstancesUNIFAC(dict):
    """Class for UNIFAC substance processing

    self['substance name'] contains :obj:`Substance`
    """
    from . import defsubs

    def __init__(self):
        pass

    def __create_phase(self,
                       subs: dict,
                       ):
        """Create phase dict for UNIFAC
        """

        # словарь веществ - групп
        tsubs = self.__get_subs(subs)
        for s in subs:
            gr = {}
            for i in tsubs[s]:
                gr[i[0]] = i[1]
            self[s] = Substance(0, gr)

    def __get_subs(self, subs: dict[str, str]) -> dict[str, list]:
        """ 'substance name': 'X*'CH2 ..' -> 'substance name': [Group, count]

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
        self.__create_phase(b)

    def get_from_dict(self, subs: dict[str, str]):
        """Load substances from input dict

        Parameters
        ----------
        subs : dict[str, str]
            Dictionary ['substance name', 'X*'CH2 ..']

        Examples
        --------
        >>> subs = {
        ...    "n-hexane": "2*CH3 4*CH2",
        ...    "butanone-2": "1*CH3 1*CH2 1*CH3CO",
        ... }
        >>> substances = SubstancesUNIFAC()
        >>> substances.get_from_dict(subs)
        """

        self.__create_phase(subs)

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
        self.__create_phase(b)
