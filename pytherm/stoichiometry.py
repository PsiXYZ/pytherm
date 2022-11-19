"""VIP module
"""
import numpy as np


def get_elements(substances) -> (np.array, np.array):
    """Extract elements from input array

    Parameters
    ----------
    substances
        Input array with substances formulas
    Returns
    -------
    np.array, np.array
        Elements array, elements matrix
    Examples
    --------
    >>> from pytherm import stoichiometry as sm
    >>> components = ['CO2', 'CO', 'H2S']
    >>> elements, elements_matrix = sm.get_elements(components)
    >>> elements
    ['C' 'O' 'H' 'S']
    >>> elements_matrix
    [[1 1 0]
     [2 1 0]
     [0 0 2]
     [0 0 1]]

    """
    elements_matrix = []  #
    elements = []  # element types list

    for s in substances:
        # поиск кол-ва элементов в системе (элемент 1 большая буква или 1б +1м)
        for i in range(len(s)):
            if s[i].isupper():
                if (i+1 < len(s)) and s[i+1].islower():
                    le = 2
                else:
                    le = 1
                if not s[i:i+le] in elements:
                    elements.append(s[i:i+le])
    # создание матрицы с колвом элементов в каждом веществе
    for i in substances:
        elements_matrix.append([0]*len(elements))
    # заполнение матрицы с колвом элементов
    for i in range(len(substances)):
        for j in range(len(substances[i])):
            # получение длины элемента (1 или 2)
            if substances[i][j].isupper():
                if (j+1 < len(substances[i])) and substances[i][j + 1].islower():
                    le = 2
                else:
                    le = 1
                # поиск цифры после элемента
                k = j + le
                str = ""
                if (k < len(substances[i])) and substances[i][k].isdigit():
                    while ((k < len(substances[i])) and substances[i][k].isdigit()):
                        str += substances[i][k]
                        k += 1
                    n = int(str)
                else:
                    n = 1
                index = elements.index(substances[i][j:j + le])
                elements_matrix[i][index] = n
    return np.array(elements), np.array(elements_matrix).T


def get_reaction_matrix(elements_matrix) -> np.array:
    r"""Generate reaction matrix from elements matrix
    Number of reaction = number of components - number of elements

    Parameters
    ----------
    elements_matrix
        Elements matrix
    Returns
    -------
    nd.array
        Reaction matrix
    Examples
    --------
    >>> from pytherm import stoichiometry as sm
    >>> components = ['CH4', 'C2H6', 'C2H4', 'C2H2']
    >>> elements, elements_matrix = sm.get_elements(components)
    >>> r_mat = sm.get_reaction_matrix(elements_matrix)
    >>> r_mat
    >>> [[ 2. -2.  1.  0.]
    ... [ 4. -3.  0.  1.]]

    """
    # необходимое кол-во реакций
    n_react = elements_matrix.shape[1] - elements_matrix.shape[0]
    # угловая матрица гаусса для el_mat
    mat = __gauss(elements_matrix)
    react_mat = np.zeros([n_react, elements_matrix.shape[1]])

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


def __gauss(m):
    """Convert input matrix to diagonal matrix
    """
    def swap(i, j):
        a = m[i]
        m[i] = m[j]
        m[j] = a

    def sort_col(c):
        for i in range(c, m.shape[0]):
            print(m[i][i])

    return m


def reaction_to_str(reaction_vector, substances, sep='*'):
    r"""Convert reaction vector to string

    Parameters
    ----------
    reaction_vector
        reaction vector
    substances
        substances array
    sep
        separator between coefficients and substances
    Returns
    -------
    str
        str reaction
    Examples
    --------
    >>> from pytherm import stoichiometry as sm
    >>> components = ['CH4', 'C2H6', 'C2H4', 'C2H2']
    >>> reaction_vector = [ 2, -2,  1,  0]
    >>> sm.reaction_to_str(reaction_vector, components)
    2*C2H6 = 2*CH4 + 1*C2H4

    """
    s1 = ""
    s2 = ""
    for i in range(len(reaction_vector)):
        if reaction_vector[i] < 0:
            s1 += str(int(abs(reaction_vector[i]))) + sep + substances[i] + " + "
        elif reaction_vector[i] > 0:
            s2 += str(int(abs(reaction_vector[i]))) + sep + substances[i] + " + "
    return s1[:-2] + "= " + s2[:-2]


def get_charge_dict():
    charge = {
        'H': +1,
        'NH4': +1,
        'Na': +1,
        'K': +1,
        'Ca': +2,
        'Mg': +2,
        'MgOH': +1,
        'HSO4': -1,
        'SO4': -2,
        'NO3': -1,
        'Cl': -1,
        'HCO3': -1,
        'CO3': -2,
        'OH': -1,
        'NH3': 0,
        'CO2 ': 0,
        'CaCO3': 0,
        'MgCO3': 0,
    }
    return charge
