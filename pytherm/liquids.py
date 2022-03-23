from math import log, log10
from .solver import minimize


def get_yinf(activity_model,
             phase: dict[str, float],
             comp_name: str) -> float:
    """Calculate activity at infinity dilution

    Returns activity at infinity dilution (y inf) for ph by extrapolation to
    low concentration

    Args:
        am (ActivityModel): model for activity calculation
        ph (dict[str, float]): _description_
        comp_name (str): name of target component from ph

    Returns:
        float: activity at infinity dilution
    """
    n_ph = 1
    n_comp = 1E-9
    x = {}
    n = {}
    for i in phase:
        n[i] = phase[i] * n_ph
    n[comp_name] = n_comp
    for i in phase:
        x[i] = n[i] / n_ph
    y = activity_model.get_y(x)
    return y[comp_name]


def get_k_molar(activity_model,
                phase1: dict[str, float],
                phase2: dict[str, float],
                comp_name: str) -> float:
    """Calculate partition coefficient using molar fraction

    Args:
        activity_model (ActivityModel): model for activity calculation
        phase1 (dict[str, float]): component dictionary for phase 1
        phase2 (dict[str, float]): component dictionary for phase 2
        comp_name (str): target component name

    Returns:
        float: partition coefficient
    """
    phase1_y = get_yinf(activity_model, phase1, comp_name)
    phase2_y = get_yinf(activity_model, phase2, comp_name)
    return phase2_y / phase1_y


def get_kp(am,
           phase1: dict[float],
           phase2: dict[float],
           molar_volume: dict[str, float],
           comp_name: str) -> float:
    """Calculate partition coefficient

    Args:
        am (ActivityModel): model for activity calculation
        phase1 (dict[float]): component dictionary for phase 1
        phase2 (dict[float]): component dictionary for phase 2
        molar_volume (dict[str, float]): molar volume dictionary
        comp_name (str): target component name

    Returns:
        float: partition coefficient
    """
    v1 = 0
    v2 = 0
    for i in phase1.keys():
        if i != comp_name:
            v1 += phase1[i] / molar_volume[i]
            v2 += phase2[i] / molar_volume[i]
    kexp = get_k_molar(am, phase1, phase2, comp_name)
    rez = kexp * v2/v1
    print("V2/V1 = ", v2/v1)
    return rez


def lle_notifier(f, ph1, ph2):
    comp = list(ph1.keys())
    n_comp = len(comp)
    # ВЫВОД
    print("F", round(log10(f)))
    n_vals = 4
    for s in comp:
        print("{} {:<12.12s} {:<7.4f} {:<7.4f}".format(
            "x", s, round(ph1[s], n_vals), round(ph2[s], n_vals)))

    s1, s2 = 0, 0
    for i in range(n_comp):
        s1 += ph1[comp[i]]
        s2 += ph2[comp[i]]
    print("{:14} {:<7.4f} {:<7.4f}".format(
        "s", round(s1, n_vals), round(s2, n_vals)))


def find_lle(phase1: dict[str, float],
             phase2: dict[str, float],
             activity_model,
             notifier=lle_notifier):
    """Calculate eqilibrium composition

    Args:
        phase1 (dict[str, float]): component dictionary for phase 1
        phase2 (dict[str, float]): component dictionary for phase 2
        activity_model (activity_model): model for activity calculation
        notifier: print results
    """
    # ph1 -> ph2
    n1 = 10  # initial amount of phase 1
    n2 = 10  # initial amount of phase 2
    min_ksi = 1e-08  # minimal ksi value
    components = list(phase1.keys())  # components list
    comp_number = len(components)  # number of components
    bounds = []  # ksi bounds for each component

    """Calculate initial amount of substance in each phase"""
    n0_1 = [0] * comp_number
    n0_2 = [0] * comp_number
    for i in range(comp_number):
        if phase1[components[i]] == 0:
            n0_1[i] = 0
            n0_2[i] = phase2[components[i]] * n2
            bounds.append([- phase2[components[i]] * n2 + min_ksi, - min_ksi])
        else:
            n0_1[i] = phase1[components[i]] * n1
            n0_2[i] = phase2[components[i]] * n2
            bounds.append([min_ksi, phase1[components[i]] * n1 - min_ksi])

    """Genereate reaction matrix for ph1 -> ph2 """
    r_mat = []
    for i in range(comp_number):
        r_mat.append([0] * comp_number)
    for i in range(comp_number):
        r_mat[i][i] = 1

    def get_loga(ksi):
        n1 = 0
        n2 = 0
        n_1 = [0] * comp_number
        n_2 = [0] * comp_number
        a1 = [0] * comp_number
        a2 = [0] * comp_number
        """calculate n for each component using ksi"""
        for i in range(comp_number):
            n_1[i] = n0_1[i] - ksi[i]
            n_2[i] = n0_2[i] + ksi[i]

        """calculate amount of phase"""
        for i in range(comp_number):
            n1 += n_1[i]
            n2 += n_2[i]

        """calculate molar fraction, х"""
        for i in range(comp_number):
            phase1[components[i]] = n_1[i] / n1
            phase2[components[i]] = n_2[i] / n2
        y1 = activity_model.get_y(phase1)
        y2 = activity_model.get_y(phase2)
        for i in range(comp_number):
            a1[i] = phase1[components[i]] * y1[components[i]]
            a2[i] = phase2[components[i]] * y2[components[i]]

        """calculate log_a for each component"""
        lg = [0] * comp_number
        for i in range(comp_number):
            lg[i] = log(a2[i] / a1[i])
        return lg

    f = minimize(get_loga, bounds)[1]
    if notifier is not None:
        lle_notifier(f, phase1, phase2)
