"""This module contain built-in datasets for UNIFAC

.. autofunction:: BIO2016_1
.. autofunction:: BIO2016_2
.. autofunction:: DOR
.. autofunction:: INF
.. autofunction:: LLE
.. autofunction:: NIST2015
.. autofunction:: PSRK
.. autofunction:: VLE
"""

from pytherm.activity.unifac import ParametersUNIFAC
from os.path import abspath

__all__ = [
    "BIO2016_1",
    "BIO2016_2",
    "DOR",
    "INF",
    "LLE",
    "NIST2015",
    "PSRK",
    "VLE",
]


def BIO2016_1() -> ParametersUNIFAC:
    """Returns BIO2016_01 :obj:`SubstancesUNIFAC` object

    Returns:
        ParametersUNIFAC: BIO2016_01 dataset
        
    References
    -----------
    Bessa2016, DOI: https://doi.org/10.1016/j.fluid.2016.05.020
    """
    params = ParametersUNIFAC(
        abspath(__file__).removesuffix("datasets.py") + "\\data\\bio2016_1.dat"
    )
    return params


def BIO2016_2() -> ParametersUNIFAC:
    """Returns BIO2016_02 :obj:`SubstancesUNIFAC` object

    Returns:
        ParametersUNIFAC: BIO2016_02 dataset
        
    References
    -----------
    Bessa2016, DOI: https://doi.org/10.1016/j.fluid.2016.05.020
    """
    params = ParametersUNIFAC(
        abspath(__file__).removesuffix("datasets.py") + "\\data\\bio2016_2.dat"
    )
    return params


def DOR() -> ParametersUNIFAC:
    """Returns DOR :obj:`SubstancesUNIFAC` object

    Returns:
        ParametersUNIFAC: DOR dataset
        
    References
    -----------
    Published DDB parameters,
    https://www.ddbst.com/PublishedParametersUNIFACDO.html
    """
    params = ParametersUNIFAC(
        abspath(__file__).removesuffix("datasets.py") + "\\data\\dor.dat"
    )
    return params


def INF() -> ParametersUNIFAC:
    """Returns INF :obj:`SubstancesUNIFAC` object

    Returns:
        ParametersUNIFAC: INF dataset
        
    References
    -----------
    Bastos1988, DOI: https://doi.org/10.1021/i200013a024
    """
    params = ParametersUNIFAC(
        abspath(__file__).removesuffix("datasets.py") + "\\data\\inf.dat"
    )
    return params


def LLE() -> ParametersUNIFAC:
    """Returns LLE :obj:`SubstancesUNIFAC` object

    Returns:
        ParametersUNIFAC: LLE dataset
        
    References
    -----------
    Magnussen1981, DOI: https://doi.org/10.1021/i200013a024
    """
    params = ParametersUNIFAC(
        abspath(__file__).removesuffix("datasets.py") + "\\data\\lle.dat"
    )
    return params


def NIST2015() -> ParametersUNIFAC:
    """Returns NIST2015 :obj:`SubstancesUNIFAC` object

    Returns:
        ParametersUNIFAC: NIST2015 dataset
        
    References
    -----------
    Kang2015, DOI: https://doi.org/10.1016/j.fluid.2014.12.042
    """
    params = ParametersUNIFAC(
        abspath(__file__).removesuffix("datasets.py") + "\\data\\nist2015.dat"
    )
    return params


def PSRK() -> ParametersUNIFAC:
    """Returns PSRK :obj:`SubstancesUNIFAC` object

    Returns:
        ParametersUNIFAC: PSRK dataset
        
    References
    -----------
    Published DDB parameters,
    https://www.ddbst.com/psrk.html
    """
    params = ParametersUNIFAC(
        abspath(__file__).removesuffix("datasets.py") + "\\data\\psrk.dat"
    )
    return params


def VLE() -> ParametersUNIFAC:
    """Returns VLE :obj:`SubstancesUNIFAC` object

    Returns:
        ParametersUNIFAC: VLE dataset
        
    References
    -----------
    Published DDB parameters,
    https://www.ddbst.com/published-parameters-unifac.html
    """
    params = ParametersUNIFAC(
        abspath(__file__).removesuffix("datasets.py") + "\\data\\vle.dat"
    )
    return params
