from pytherm.activity.sit import SIT
from pytherm.systems_test.chemicalsystem import ChemicalReaction, SolidPhase
from pytherm.systems_test.equilibriumsystem import ElectrolyteSolution, EquilibriumSystem, EquilibriumSolver, SolubilitySP

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

    ph_solid = [
        'Cocl2s',
    ]
    reactions_sol = {
        ChemicalReaction(
            log10_k=1,
            reaction_str="1*CoCl2s = 1*Co_+2 + 2*Cl_-1"
        ),
    }

    system = EquilibriumSystem()

    solution = ElectrolyteSolution(substances_dict=ph, reactions_list=reactions_list)
    activity_model = SIT(solution.substances_list)
    solution.set_activity_model(activity_model)
    system.add_phase(solution)

    # solid = SolidPhase(ph_solid)

    # solubility_exchange = SolubilitySP(solid, solution)
    # solubility_exchange.set_reactions(reactions_sol)

    solver = EquilibriumSolver()
    system.set_solver(solver)



    system.equilibrate(298)

    conc = solution.get_concentrations(solver.ksi)
    subs = solution.substances_list
    for i in range(len(conc)):
        print(f"{subs[i]}, {conc[i]:.6f}")


# import timeit
# n=50
# elapsed_time = timeit.timeit(test, number=n)/n
# print('Elapsed time: ', elapsed_time)

test()