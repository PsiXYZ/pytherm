import os

path = os.path.abspath(__file__)[:-5:]

class Sub:
    def __init__(self, conc, groups):
        if conc != None:
            self.x = float(conc)
        self.groups = groups  # словарь {Имя группы : кол-во}


class Group:
    def __init__(self, id, r, q):
        self.id = int(id)
        self.R = float(r)
        self.Q = float(q)

# return a(i,j) matrix
def get_inter(mode):
    df = open(path + mode + "\\2t.csv")
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
        table[int(i[0]) - 1][int(i[1]) - 1] = [float(i[2]), float(i[3]), float(i[4])]
    for i in range(n):
        table[i][i] = [0, 0, 0]
    return table

# return {gr name: Group}
def get_groups(mode):
    df = open(path + mode + "\\1t.csv")
    buf = []
    for line in df:
        buf.append(line.split(";"))
    df.close()

    data = {}
    for i in range(1, len(buf)):
        data[buf[i][1]] = Group(buf[i][0], buf[i][2], buf[i][3])
        # data[buf[i][1]] = [buf[i][2], buf[i][3], buf[i][0]]

    return data


# возвращает словарь  {название вещества: [группа, кол-во]}, вкачивает файл полностью, лучше потом сделать поиск
def get_tsubs(mode, subs_mode):

    data = {}
    if subs_mode == "с":
        df = open(path + "\\sub.csv")
    else:
        df = open(path + "\\" + mode + "\\sub.csv")

    buf = []
    for line in df:
        buf.append(line.split(":"))
    for i in range(len(buf)):
        buf[i][1] = buf[i][1].replace("\n", "")
        buf[i][1] = buf[i][1][1::]

    for i in buf:
        print(i)
        b = i[1].split(" ")
        a = []
        for j in b:
            a.append([j.split("*")[1], int(j.split("*")[0]) ])
        data[i[0]] = a

    return data


# {name: x} --> {name: sub}
def create_ph(inp, mode, subs_mode):
    phase = {}
    tsubs = get_tsubs(mode, subs_mode)  # словарь веществ - групп

    for s in inp:
        gr = {}
        for i in tsubs[s]:
            gr[i[0]] = i[1]
        phase[s] = Sub(0, gr)

    return phase