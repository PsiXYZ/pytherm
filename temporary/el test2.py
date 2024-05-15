from pytherm.activity.sit import SIT
from pytherm.systems.chemicalsystem import ChemicalReaction
from pytherm.systems.phase import ElectrolyteSolution
from pytherm.systems.equilibriumsystem import EquilibriumSystem

def test():
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

    system = EquilibriumSystem()

    solution = ElectrolyteSolution()
    solution.add_reactions(reactions_list)
    solution.add_substances(list(ph.keys()))

    activity_model = SIT(solution.substances_str_list)
    solution.set_activity_model(activity_model)

    solution.set_concentrations(ph)

    system.add_phase(solution)

    system.equilibrate(298)

    conc = solution.conc
    subs = solution.substances_str_list
    for i in range(len(conc)):
        print(f"{subs[i]}, {conc[i]:.6f}")

# import timeit
# n=50
# elapsed_time = timeit.timeit(test, number=n)/n
# print('Elapsed time: ', elapsed_time)

test()