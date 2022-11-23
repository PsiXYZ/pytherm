from .parametersIO import ParametersPitzerNew
import numpy as np


# pitzer.dat
f = lambda T: np.array([1,
                        1/T - 1/298.15,
                        np.log(T/298.15),
                        T - 298.15,
                        T**2 - 298.15 ** 2,
                        1/(T**2) - 1/(298.15**2)])
pitzer_dataset = ParametersPitzerNew()
pitzer_dataset.load_from_dat('./pitzer.dat', f)
