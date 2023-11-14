from pytherm.activity.unifac import ParametersUNIFAC
from os.path import abspath

__all__ = [
    "get_BIO2016_1",
    "get_BIO2016_2",
    "get_DOR",
    "get_INF",
    "get_LLE",
    "get_NIST2015",
    "get_PSRK",
    "get_VLE",
]


def get_BIO2016_1() -> ParametersUNIFAC:
    params = ParametersUNIFAC(
        abspath(__file__).removesuffix("datasets.py") + "\\data\\bio2016_1.dat"
    )
    return params


def get_BIO2016_2() -> ParametersUNIFAC:
    params = ParametersUNIFAC(
        abspath(__file__).removesuffix("datasets.py") + "\\data\\bio2016_2.dat"
    )
    return params


def get_DOR() -> ParametersUNIFAC:
    params = ParametersUNIFAC(
        abspath(__file__).removesuffix("datasets.py") + "\\data\\dor.dat"
    )
    return params


def get_INF() -> ParametersUNIFAC:
    params = ParametersUNIFAC(
        abspath(__file__).removesuffix("datasets.py") + "\\data\\inf.dat"
    )
    return params


def get_LLE() -> ParametersUNIFAC:
    params = ParametersUNIFAC(
        abspath(__file__).removesuffix("datasets.py") + "\\data\\lle.dat"
    )
    return params


def get_NIST2015() -> ParametersUNIFAC:
    params = ParametersUNIFAC(
        abspath(__file__).removesuffix("datasets.py") + "\\data\\nist2015.dat"
    )
    return params


def get_PSRK() -> ParametersUNIFAC:
    params = ParametersUNIFAC(
        abspath(__file__).removesuffix("datasets.py") + "\\data\\psrk.dat"
    )
    return params


def get_VLE() -> ParametersUNIFAC:
    params = ParametersUNIFAC(
        abspath(__file__).removesuffix("datasets.py") + "\\data\\vle.dat"
    )
    return params
