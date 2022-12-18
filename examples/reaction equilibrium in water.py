from pytherm.systems import *
from pytherm.activity.sit import SIT

ph = {
    'Co_+2': 0.1,
    'Cl_-1': 0.7,
    'Na_+1': 0.5
}

reactions_list = [
    ChemicalReaction(
        log10_k=0.570,
        reaction_str="1*Co_+2 + 1*Cl_-1 = 1*CoCl_+1"
    ),
    ChemicalReaction(
        log10_k=0.020,
        reaction_str="1*Co_+2 + 2*Cl_-1 = 1*CoCl2"
    ),
    ChemicalReaction(
        log10_k=-1.710,
        reaction_str="1*Co_+2 + 3*Cl_-1 = 1*CoCl3_-1"
    ),
    ChemicalReaction(
        log10_k=-2.090,
        reaction_str="1*Co_+2 + 4*Cl_-1 = 1*CoCl4_-2"
    ),
]

chem_sys = EquilibriumSystem(
    list(ph.keys()),
    reactions_list,
)
chem_sys.set_concentrations(ph)
chem_sys.set_activity_model(SIT(chem_sys.substances))

solver = EquilibriumSolver(chem_sys)
solver.equilibrate()

conc = chem_sys.get_concentrations(solver.ksi)
subs = chem_sys.substances
for i in range(len(conc)):
    print(f"{subs[i]}, {conc[i]:.6f}")
