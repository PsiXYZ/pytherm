import os

path = os.path.abspath(__file__)[:-5:]


class Substance:
    """Class for unifac substance using Group
    """

    def __init__(self, concentration, groups):
        if concentration != None:
            self.x = float(concentration)
        self.groups = groups  # словарь {Имя группы : кол-во}


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


def get_interactions(unifac_mode: str) -> list:
    """Gets interaction parameters from unifac_mode file

    Args:
        unifac_mode (str): unifac parameters (VLE, LLE etc.)

    Returns:
        list: a(i, j) interaction matrix
    """
    df = open(path + unifac_mode + "\\2t.csv")
    n = int(df.readline()[:-1].split(";")[0])

    table = []
    for i in range(n):
        a = []
        for j in range(n):
            a.append("-")
        table.append(a)

    buf = []
    for line in df:
        buf.append(line[:-1:].split(";"))

    for i in buf[1::]:
        table[int(i[0]) - 1][int(i[1]) - 1] = [float(i[2]),
                                               float(i[3]), float(i[4])]
    for i in range(n):
        table[i][i] = [0, 0, 0]
    return table


def get_groups(unifac_mode: str) -> dict[str, Group]:
    """Create dictionary with Group from unifac_mode file

    Args:
        unifac_mode (str): unifac parameters(VLE, LLE etc.)

    Returns:
        dict: {gr name: Group}
    """
    df = open(path + unifac_mode + "\\1t.csv")
    buf = []
    for line in df:
        buf.append(line.split(";"))
    df.close()

    data = {}
    for i in range(1, len(buf)):
        data[buf[i][1]] = Group(buf[i][0], buf[i][2], buf[i][3])
    return data


def get_tsubs(unifac_mode: str, substance_source: str) -> dict:
    """

    Args:
        unifac_mode (str): unifac parameters (VLE, LLE etc.)
        substance_source (str): source for substance group
            ('general' for sub from db, 'unifac' for for sub from mode\db)

    Returns:
        dict: {substance name: [Group, count]}
    """
    data = {}
    if substance_source == "general":
        df = open(path + "\\sub.csv")
    else:
        df = open(path + "\\" + unifac_mode + "\\sub.csv")

    buf = []
    for line in df:
        buf.append(line.split(":"))
    for i in range(len(buf)):
        buf[i][1] = buf[i][1].replace("\n", "")
        buf[i][1] = buf[i][1][1::]

    for i in buf:
        b = i[1].split(" ")
        a = []
        for j in b:
            a.append([j.split("*")[1], float(j.split("*")[0])]) #int -> float
        data[i[0]] = a

    return data


def create_phase(inp: dict,
                 unifac_mode: str,
                 substance_source: str) -> dict[str, Substance]:
    """Create phase dict for Unifac

    Args:
        inp (dict): {substance_name: concentration}  dictionary
        unifac_mode (str): unifac parameters(VLE, LLE etc.)
        substance_source (str): source for substance group
            ('general' for sub from db, 'unifac' for for sub from mode\db)

    Returns:
        dict[str, Substance]: dict for Unifac calculation
    """
    phase = {}
    tsubs = get_tsubs(unifac_mode, substance_source)  # словарь веществ - групп

    for s in inp:
        gr = {}
        for i in tsubs[s]:
            gr[i[0]] = i[1]
        phase[s] = Substance(0, gr)
    return phase
