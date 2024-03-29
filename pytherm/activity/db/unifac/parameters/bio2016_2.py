from .. import ParametersUNIFAC

t1 = [
    [1, 'CH3', 0.9011, 0.848],
    [1, 'CH2', 0.6744, 0.54],
    [1, 'CH', 0.4469, 0.228],
    [1, 'C', 0.2195, 0.0],
    [2, 'CH=CH', 1.1167, 0.867],
    [3, 'OH', 1.0, 1.2],
    [4, 'H2O', 0.92, 1.4],
    [5, 'CH2COO', 1.6764, 1.42],
    [6, 'OHgly', 1.0, 1.2],
]

t2 = [
    [1, 1, 0.0, 0, 0],
    [1, 2, 380.42, 0, 0],
    [1, 3, 394.21, 0, 0],
    [1, 4, 2577.21, 0, 0],
    [1, 5, 1127.95, 0, 0],
    [1, 6, 1170.1, 0, 0],
    [2, 1, 2037.3, 0, 0],
    [2, 2, 0.0, 0, 0],
    [2, 3, 1192.65, 0, 0],
    [2, 4, 74.65, 0, 0],
    [2, 5, -12.05, 0, 0],
    [2, 6, 91.75, 0, 0],
    [3, 1, 134.05, 0, 0],
    [3, 2, 266.61, 0, 0],
    [3, 3, 0.0, 0, 0],
    [3, 4, 14.68, 0, 0],
    [3, 5, 202.16, 0, 0],
    [3, 6, 0.0, 0, 0],
    [4, 1, 129.97, 0, 0],
    [4, 2, 351.07, 0, 0],
    [4, 3, -78.61, 0, 0],
    [4, 4, 0.0, 0, 0],
    [4, 5, -7.2, 0, 0],
    [4, 6, 85.32, 0, 0],
    [5, 1, -446.43, 0, 0],
    [5, 2, 4.763, 0, 0],
    [5, 3, 555.53, 0, 0],
    [5, 4, 326.49, 0, 0],
    [5, 5, 0.0, 0, 0],
    [5, 6, -21.65, 0, 0],
    [6, 1, 375.35, 0, 0],
    [6, 2, 1720.1, 0, 0],
    [6, 3, 0.0, 0, 0],
    [6, 4, -8.415, 0, 0],
    [6, 5, -23.25, 0, 0],
    [6, 6, 0.0, 0, 0],
]
BIO2016_2 = ParametersUNIFAC("BIO2016_2")
BIO2016_2.set_type("classic")
BIO2016_2.set_comb(t1)
BIO2016_2.set_res(t2)
