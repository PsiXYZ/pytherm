from pytherm.systems.chemicalsystem import EquilibriumReaction, EquilibriumSystem
import numpy as np

def create_system():
    ph1 = {
        'A{l1}': 1.0,
        'B{l1}': 0.0,
        'C{l1}': 0.0,
        'D{l1}': 0.0,
        'E{l1}': 0.0,
        'F{l1}': 0.0,
        'G{l1}': 0.0,
    }
    ph2 = {
        'A{l2}': 0.0,
        'B{l2}': 0.0,
        'C{l2}': 1.0,
        'D{l2}': 0.0,
        'E{l2}': 0.0,
        'F{l2}': 0.0,
        'G{l2}': 0.0,
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
        EquilibriumReaction(
            log10_k=np.log10(0.1),
            reaction_str="1*A{l1} = 1*D{l1}"
        ),
        EquilibriumReaction(
            log10_k=np.log10(0.1),
            reaction_str="1*A{l1} = 1*E{l1}"
        ),
        EquilibriumReaction(
            log10_k=np.log10(0.01),
            reaction_str="1*A{l1} = 1*F{l1}"
        ),
        EquilibriumReaction(
            log10_k=np.log10(0.01),
            reaction_str="1*A{l1} = 1*G{l1}"
        ),
        # EquilibriumReaction(
        #     log10_k=np.log10(0.01),
        #     reaction_str="1*D{l1} = 1*E{l1}"
        # ),
        EquilibriumReaction(
            log10_k=np.log10(1.0),
            reaction_str="1*B{l1} = 1*B{l2}"
        ),
        EquilibriumReaction(
            log10_k=np.log10(1.0),
            reaction_str="1*D{l1} = 1*D{l2}"
        ),
    ]
    # print(*reactions_list)
    system = EquilibriumSystem()
    # system.init_by_reactions(reactions_list)
    system.add_phase_by_dict(ph1)
    system.add_phase_by_dict(ph2)
    system.add_reactions_by_list(reactions_list)
    system.assemble()

    return system, ph

def test0():
    system, ph = create_system()
    system.equilibrate(ph, solver_type="old")
    print("old")
    print(system.ksi)
    print(system.get_Q(system.ksi))

def test1():
    system, ph = create_system()
    system.equilibrate(ph, solver_type="scipy")
    print("scipy")
    print(system.ksi)
    print(system.get_Q(system.ksi))

def test2():
    system, ph = create_system()
    system.equilibrate(ph, solver_type="custom_lm")
    print("custom")
    print(system.ksi)
    print(system.get_Q(system.ksi))  

import timeit
n=100
elapsed_time = timeit.timeit(test2, number=n)/n
print('Elapsed time: ', elapsed_time)

# test0()
# test1()
# test2()