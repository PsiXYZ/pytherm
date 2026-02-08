"""
This module contains classes and methods to calculate activity coefficients
with the UNIQUAC model.

How to use
----------
To use UNIQUAC it's nessesary to set-up the model:
    >>> from pytherm.activity import uniquac as uq
    >>> T = 25.0 + 273.15
    >>> xs = [0.7273, 0.0909, 0.1818]
    >>> rs = [0.92, 2.1055, 3.1878]
    >>> qs = [1.4, 1.972, 2.4]
    >>> inter = [
    ...     [[0, 0], [0, 526.02], [0, 309.64]],
    ...     [[0, -318.06], [0, 0], [0, -91.532]],
    ...     [[0, 1325.1], [0, 302.57], [0, 0]],
    ... ]
    >>> am = uq.UNIQUAC(rs, qs, inter)
    >>> print(am.get_y(xs, T))
    [1.570393443107605, 0.29482409358024597, 18.11433219909668]


UNIQUAC class
------------
.. autoclass:: UNIQUAC
    :members: __init__, get_a, get_y
   
"""

from __future__ import annotations

from pytherm.cpp import (
    ActivityModel,
    UNIQUAC,
)

__all__ = [
    "ActivityModel",
    "UNIQUAC",
]
