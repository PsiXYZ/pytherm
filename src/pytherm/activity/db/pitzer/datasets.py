from .parametersIO import ParametersPitzerNew
import numpy as np
import os

path = os.path.dirname(os.path.realpath(__file__))

# pitzer.dat
f = lambda T: np.array([1,
                        1/T - 1/298.15,
                        np.log(T/298.15),
                        T - 298.15,
                        T**2 - 298.15 ** 2,
                        1/(T**2) - 1/(298.15**2)
                        ])
pitzer_dataset = ParametersPitzerNew()
pitzer_dataset.load_from_dat(path + '/pitzer.dat', f)
