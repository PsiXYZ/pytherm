from pytherm import base, stoichiometry, gaseq
from pytherm import gaseq as gaseq2


T = 900 + 273  # [K]
P = 1
system = {
    # 'CO2': 2,
    # 'CO': 0.1,
    # 'H2O': 6,
    'H2': 0.0000,
    'CH4': 4,
    'C2H6': 1,
    'C2H4': 0.0,
    'C2H2': 0.0,
    'C6H6': 0.0,
}
hf = {
    'H2': 0,
    'CO': -110.53,
    'CO2': -393.51,
    'H2O': -241.81,
    'CH4': -74.85,
    'C2H6': -84.67,
    'C2H4': 52.3,
    'C2H2': 226.75,
    'C3H8': -103.85,
    'C3H6': 20.41,
    'C6H6': 82.93,
    'CH4O': -201,
    'C2H6O': -234.8,
}
sf = {
    'H2': 130.52,
    'CO': 197.55,
    'CO2': 213.66,
    'H2O': 188.72,
    'CH4': 186.27,
    'C2H6': 229.49,
    'C2H4': 219.45,
    'C2H2': 200.82,
    'C3H8': 269.91,
    'C3H6': 266.94,
    'C6H6': 269.2,
    'CH4O': 239.76,
    'C2H6O': 281.38,

}
kelly = {
    'H2': (27.28, 3.26, 0, 0.5),
    'CO': (28.41, 4.1, 0, -0.46),
    'CO2': (44.14, 9.04, 0, -8.54),
    'H2O': (30, 10.71, 0, 0.33),
    'CH4': (14.32, 74.66, -17.43, 0),
    'C2H6': (5.75, 175.11, -57.85, 0),
    'C2H4': (11.32, 122.01, -37.9, 0),
    'C2H2': (20.44, 66.65, -26.48, 0),
    'C3H8': (1.72, 270.75, -94.48, 0),
    'C3H6': (12.44, 188.38, -47.6, 0),
    'C6H6': (-21.09, 400.12, -169.87, 0),
    'CH4O': (15.28, 105.2, -31.04, 0),
    'C2H6O': (10.99, 204.7, -74.2, 0),
}
components = tuple(system.keys())

el, e_su = sm.get_elements(components)
r_mat = sm.get_reaction_matrix(e_su)

print("components:", components)
print("elements:",  el)
print("elements matrix:\n", e_su)
print(f"reaction matrix:\n{components}\n{r_mat}")

G = {}
for component in components:
    G[component] = base.get_Gt_kelly(
        h0=hf[component],
        s0=sf[component],
        kelly=kelly[component],
        T=T
    )
print("G", G)

K = base.get_k(
    G=tuple(G.values()),
    T=T,
    r_mat=r_mat,
)
print("K", K)

for i in r_mat:
    print(sm.reaction_to_str(i, components))

out = gaseq2.find_eq(
    n0=tuple(system.values()),
    r_mat=r_mat,
    K=K,
    T=T,
    P=P
)

n0 = 0
for i in system:
    n0 += system[i]
for i in system:
    system[i] /= n0
print(system)

for i in range(len(components)):
    system[components[i]] = out[i]
n0 = 0
for i in system:
    n0 += system[i]
for i in system:
    system[i] /= n0
print(system)