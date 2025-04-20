=======
Pytherm
=======

.. image:: https://readthedocs.org/projects/pytherm/badge/?version=latest
    :target: https://pytherm.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. contents::

What is Pytherm?
----------------
Pytherm is an open-source scientific Python package for thermodynamic modeling. 

Library include ready to use functions for VLE and LLE modeling.
Modeling can be provided by different built-in EOS and activity models with 
user-friendly interfaces.

Current features
-----------------
* Activity models:
 
  * UNIFAC and Modified UNIFAC (6 parameters sets)
  * UNIQUAC
  * PITZER (full form)
  * Specific ion Interaction

* EOS models:

  * PSRK 

* Calculations:

  * Solubility in 2 and 3 component systems 
  * LLE for non-electrolyte systems 
  * VLE using EOS
  * Multicomponent equilibrium calculations using PITZER ans SIT models

Getting Started
---------------
Example of activity coefficient calculations using UNIFAC:

.. code-block:: python

    import pytherm.activity.unifac as uf 
    subs = {
        "n-hexane": "2*CH3 4*CH2", 
        "butanone-2": "1*CH3 1*CH2 1*CH3CO",
    }
    system = {
        'n-hexane': 0.5,
        'butanone-2': 0.5,
    }

    substances = uf.datasets.SubstancesUNIFAC()
    substances.get_from_dict(subs)

    am = uf.Unifac(dataset=uf.datasets.DOR, substances=substances)
    print(am.get_y(system=system, T=298))

Some more examples of use can be found in /examples

Installation
-------------
  $ git clone git://github.com/PsiXYZ/pytherm
  $ git submodule update --init --recursive
  $ pip install .

Documentation
--------------
Unfortunately at the moment the documentation is broken.
    https://pytherm.readthedocs.io/en/latest/index.html

Roadmap
-------
Future features:

* Activity models:

  * COSMO-SAC
  * NRTL
  * Regular solution
  * LIFAC
  * LIQUAC

* Calculations:

  * Convex hull analysis for phase equilibrium
