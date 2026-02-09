from pytherm.systems.chemicalsystem import EquilibriumReaction, EquilibriumSystem
import numpy as np


def test():
    ph1 = {
        'A{l1}': 1.0,
        'B{l1}': 0.0,
        'C{l1}': 0.0
    }
    ph2 = {
        'A{l2}': 0.0,
        'B{l2}': 0.0,
        'C{l2}': 1.0
    }

    ph = {
        'l1': ph1,
        'l2': ph2
    }

    reactions_list = [
        EquilibriumReaction(
            log10_k=np.log10(0.1),
            reaction_str="1*A{l1} = 1*B{l1}"
        ),
        # EquilibriumReaction(
        #     log10_k=0.0,
        #     reaction_str="1*A{l1} = 1*A{l2}"
        # ),
        EquilibriumReaction(
            log10_k=np.log10(1.0),
            reaction_str="1*B{l1} = 1*B{l2}"
        ),
    ]
    # print(*reactions_list)

    system = EquilibriumSystem()
    # system.init_by_reactions(reactions_list)
    system.add_phase_by_dict(ph1)
    system.add_phase_by_dict(ph2)
    system.add_reactions_by_list(reactions_list)
    system.assemble()
    system.equilibrate(ph)

    
    # ksi = np.array([0.1, 0.2, 0.3])
    # system.get_p(ksi)
    
    


    # system = EquilibriumSystem()

    # solution = ElectrolyteSolution()
    # solution.add_reactions(reactions_list)
    # solution.add_substances(list(ph.keys()))

    # activity_model = SIT(solution.substances_str_list)
    # solution.set_activity_model(activity_model)

    # solution.set_concentrations(ph)

    # system.add_phase(solution)

    # system.equilibrate(298)

    # conc = solution.conc
    # subs = solution.substances_str_list
    # for i in range(len(conc)):
    #     print(f"{subs[i]}, {conc[i]:.6f}")

# import timeit
# n=50
# elapsed_time = timeit.timeit(test, number=n)/n
# print('Elapsed time: ', elapsed_time)

test()