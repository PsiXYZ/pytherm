"""
This module contains classes and methods to calculate activity coefficients
with the UNIFAC model.

How to use
----------
To use UNIFAC it's nessesary to define substances with :obj:`.SubstancesUNIFAC`,
load parameters with :obj:`.ParametersUNIFAC` and set-up the model:
    >>> import pytherm.activity.unifac as uf
    >>> subs = {
    ...    "n-hexane": "2*CH3 4*CH2",
    ...    "butanone-2": "1*CH3 1*CH2 1*CH3CO",
    ... }
    >>> system = {
    ...    'n-hexane': 0.5,
    ...    'butanone-2': 0.5,
    ... }
    >>> substances = uf.datasets.SubstancesUNIFAC()
    >>> substances.get_from_dict(subs)
    >>> am = uf.UNIFAC(dataset=uf.datasets.DOR, substances=substances)
    >>> am.get_y(conc=system, T=298)
    {'n-hexane': 1.514766775270851, 'butanone-2': 1.4331647782163541}

UNIFAC class
------------
.. autoclass:: UNIFAC
    :members: __init__, get_a, get_y
   
UNIFAC_W class
---------------
.. autoclass:: UNIFAC_W
    :members: __init__, get_a, get_y
    :show-inheritance:
    :member-order: bysource

SubstancesUNIFAC class
----------------------
Substances must be a special :obj:`.SubstancesUNIFAC` object

.. autoclass:: SubstancesUNIFAC
    :members: get_from_dict

ParametersUNIFAC
-----------------
Parameters must be a special :obj:`.ParametersUNIFAC` object

.. autoclass:: ParametersUNIFAC
    :members: __init__

Build-in parameters
-------------------
There are some ready to use :obj:`.ParametersUNIFAC` objects in :obj:`.unifac.datasets`:

* Classic UNIFAC:
    * :obj:`.unifac.datasets.VLE` (:obj:`pytherm.parameters.unifac.datasets.VLE`) [1]_
    * :obj:`.unifac.datasets.LLE` (:obj:`pytherm.parameters.unifac.datasets.LLE`) [2]_
    * :obj:`.unifac.datasets.INF` (:obj:`pytherm.parameters.unifac.datasets.INF`) [3]_
    * :obj:`.unifac.datasets.BIO2016_1` (:obj:`pytherm.parameters.unifac.datasets.BIO2016_1`) [6]_
    * :obj:`.unifac.datasets.BIO2016_2` (:obj:`pytherm.parameters.unifac.datasets.BIO2016_2`) [6]_
* Modified UNIFAC:
    * :obj:`.unifac.datasets.DOR` (:obj:`pytherm.parameters.unifac.datasets.DOR`) [4]_
    * :obj:`.unifac.datasets.NIST2015` (:obj:`pytherm.parameters.unifac.datasets.NIST2015`) [5]_

References
-----------

.. [1] Published DDB parameters, 2021 JAN,
    https://www.ddbst.com/published-parameters-unifac.html
.. [2] Magnussen1981, DOI: https://doi.org/10.1021/i200013a024
.. [3] Bastos1988, DOI: https://doi.org/10.1021/i200013a024
.. [4] Published DDB parameters, 2021 JAN,
    https://www.ddbst.com/PublishedParametersUNIFACDO.html
.. [5] Kang2015, DOI: https://doi.org/10.1016/j.fluid.2014.12.042
.. [6] Bessa2016, DOI: https://doi.org/10.1016/j.fluid.2016.05.020
"""

from __future__ import annotations


from pytherm.cpp import (
    ParametersUNIFAC,
    SubstancesUNIFAC,
    UNIFAC,
    UNIFAC_W,
    ActivityModel,
)

import pytherm.parameters.unifac as datasets

__all__ = [
    "ParametersUNIFAC",
    "SubstancesUNIFAC",
    "UNIFAC",
    "UNIFAC_W",
    "ActivityModel",
    "datasets",
]
