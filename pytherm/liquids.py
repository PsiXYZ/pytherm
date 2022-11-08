from math import log, log10
from pytherm import constants
from .solver import minimize
from pytherm.activity.activity_model import ActivityModel

R = constants.R


def get_yinf(activity_model: ActivityModel,
             system: dict[str, float],
             comp_name: str) -> float:
    """Calculate activity coefficient at infinity dilution

    Parameters
    ----------
    activity_model : ActivityModel
        model for activity calculation
    phase : dict[str, float]
        Input dictionary {"Substance name": concentration}
    comp_name : str
        name of target component from ph

    Returns
    -------
    float
        activity coefficient at infinity dilution
    """
    n_ph = 1
    n_comp = 1E-9
    x = {}
    n = {}
    for i in system:
        n[i] = system[i] * n_ph
    n[comp_name] = n_comp
    for i in system:
        x[i] = n[i] / n_ph
    y = activity_model.get_y(x)
    return y[comp_name]


def get_k_molar(activity_model: ActivityModel,
                system1: dict[str, float],
                system2: dict[str, float],
                comp_name: str) -> float:
    """Calculate partition coefficient using molar fraction

    Parameters
    ----------
    activity_model : ActivityModel
        model for activity calculation
    system1 : dict[str, float]
        component dictionary for phase 1
    phase2 : dict[str, float]
        component dictionary for phase 2
    comp_name : str
        target component name

    Returns
    -------
    float
        partition coefficient
    """
    phase1_y = get_yinf(activity_model, system1, comp_name)
    phase2_y = get_yinf(activity_model, system2, comp_name)
    return phase2_y / phase1_y


def get_kp(activity_model: ActivityModel,
           system1: dict[str, float],
           system2: dict[str, float],
           molar_volumes: dict[str, float],
           comp_name: str) -> float:
    """Calculate partition coefficient

    Parameters
    ----------
    activity_model : ActivityModel
        model for activity calculation
    system1 : dict[str, float]
        component dictionary for phase 1
    system2 : dict[str, float]
        component dictionary for phase 2
    molar_volumes : dict[str, float]
        molar volume dictionary
    comp_name : str
        target component name

    Returns
    -------
    float
        partition coefficient
    """
    v1 = 0
    v2 = 0
    for i in system1.keys():
        if i != comp_name:
            v1 += system1[i] / molar_volumes[i]
            v2 += system2[i] / molar_volumes[i]
    kexp = get_k_molar(activity_model, system1, system2, comp_name)
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
             temperature=298,
             notifier=lle_notifier):
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
        y1 = activity_model.get_y(phase1, temperature=temperature)
        y2 = activity_model.get_y(phase2, temperature=temperature)
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


def calculate_solubility(activity_model: ActivityModel,
                         comp_name,
                         T: float,
                         T_m: float,
                         H_m: float,
                         bounds=(1e-20, 0.99),
                         ftol=1e-18,
                         fabs=1e-10):
    phase = {}
    phase[comp_name[0]] = 1
    phase[comp_name[1]] = 0
    lnx = H_m / R * (1/T_m - 1/T)

    a = bounds[0]
    b = bounds[1]

    def f(x):
        phase[comp_name[1]] = x
        phase[comp_name[0]] = 1 - x
        y = activity_model.get_y(phase, T=T)
        return (log(phase[comp_name[1]] * y[comp_name[1]]) - lnx)

    while(1):
        x = (a + b) / 2
        if f(x) * f(a) < 0:
            b = x
        elif f(x) * f(b) < 0:
            a = x
        if abs(b - a) < ftol or abs(f(x)) < fabs:
            break
    return x
