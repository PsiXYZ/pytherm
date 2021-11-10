import numpy as np


def get_el(su):
    """
    По листу веществ возвращает два листа - первый это лист с набором элементов, второй это двумерный список с колвом
    элементов в каждом веществе

    :param su: substance list
    :return: 2 list - element list (C,H,O..), number of elements in sub ()
    """

    e_su = []  #
    el = []  # element types list

    for s in su:
        # поиск кол-ва элементов в системе (элемент 1 большая буква или 1б +1м)
        for i in range(len(s)):
            if s[i].isupper():
                if (i+1 < len(s)) and s[i+1].islower():
                    le = 2
                else:
                    le = 1
                if not s[i:i+le] in el:
                    el.append(s[i:i+le])
    # создание матрицы с колвом элементов в каждом веществе
    for i in su:
        e_su.append([0]*len(el))
    # заполнение матрицы с колвом элементов
    for i in range(len(su)):
        for j in range(len(su[i])):
            # получение длины элемента (1 или 2)
            if su[i][j].isupper():
                if (j+1 < len(su[i])) and su[i][j+1].islower():
                    le = 2
                else:
                    le = 1
                # поиск цифры после элемента
                k = j + le
                str = ""
                if (k < len(su[i])) and su[i][k].isdigit():
                    while ((k < len(su[i])) and su[i][k].isdigit()):
                        str += su[i][k]
                        k += 1
                    n = int(str)
                else:
                    n = 1
                index = el.index(su[i][j:j+le])
                e_su[i][index] = n
    return [el, np.array(e_su).T]


# по матрице элементов составляет реакции
def get_reaction(el_mat):
    # необходимое кол-во реакций
    n_react = el_mat.shape[1] - el_mat.shape[0]
    # угловая матрица гаусса для el_mat
    mat = gauss(el_mat)
    react_mat = np.zeros([n_react, el_mat.shape[1]])

    print(n_react)
    print(el_mat)

    # поиск колонки с детерминантом не равным 0
    i_col = 0
    for i in range(mat.shape[1]):
        det = np.linalg.det(mat[0:mat.shape[0], i:mat.shape[0] + i])
        if det != 0:
            i_col = i
            break

    r = 0
    a = mat[0:mat.shape[0], i_col:mat.shape[0] + i_col]
    for i in range(mat.shape[1]):
        if i not in range(i_col, i_col + mat.shape[0]):
            b = np.zeros([mat.shape[0], 1])
            b = - mat[:, i]
            x = np.linalg.solve(a, b)

            for j in range(x.shape[0]):
                react_mat[r, j + i_col] = x[j]
            react_mat[r, i] = 1
            r += 1
    return react_mat


def gauss(m):
    def swap(i, j):
        a = m[i]
        m[i] = m[j]
        m[j] = a

    def sort_col(c):
        for i in range(c, m.shape[0]):
            print(m[i][i])

    return m


# строковый вывод el
def el_string(su, el, e_su):
    print("___el test___")
    for i in range(len(su)):
        s = ""
        s += su[i] + " = "
        for j in range(len(el)):
            s += str(e_su[i][j]) + str(el[j]) + " + "
        print(s[0:-2:])
    print("_____________")


def react_string(c_v, su):
    s1 = ""
    s2 = ""
    for i in range(len(c_v)):
        if c_v[i] < 0:
            s1 += str(int(abs(c_v[i]))) + "•" + su[i] + " + "
        elif c_v[i] > 0:
            s2 += str(int(abs(c_v[i]))) + "•" + su[i] + " + "
    return s1[:-2] + "= " + s2[:-2]