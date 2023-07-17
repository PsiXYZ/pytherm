from pytherm.eos import PSRK
import pytherm.activity.unifac as uf
from pytherm import vle

ms_params = {
    'CH4': {
        'c1': 0.4926,
        'c2': 0,
        'c3': 0
    },
    'CO2': {
        'c1': 0.8252,
        'c2': 0.2515,
        'c3': -1.7039
    },
}
cr_params = {
    'CH4': {
        'Tc': 190.6,
        'Pc': 46 * 100000,
    },
    'CO2': {
        'Tc': 304.2,
        'Pc': 73.8 * 100000,
    },
}
subs = {
    "CO2": "1*CO2",
    "CH4": "1*CH4",
}

T = 219
xi = 0.7
system = {
    'CH4': xi,
    'CO2': 1-xi,
}

substances = uf.datasets.SubstancesUNIFAC()
substances.get_from_dict(subs)

am = uf.UNIFAC(dataset=uf.datasets.PSRK,
               substances=substances,
               dict_mode=True)
model = PSRK(system=system, cr_params=cr_params, 
             ms_params=ms_params, activity_model=am)

bubble_curve, dew_curve = vle.fit_Pxy(model=model, T=T)



import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.plot(bubble_curve[0], bubble_curve[1])
ax.plot(dew_curve[0], dew_curve[1])

# ax.set_ylim([0, 75])
ax.set(xlabel='x, y', ylabel='P (Pa)')
ax.grid()

plt.show()
