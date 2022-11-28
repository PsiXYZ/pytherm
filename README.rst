=======
Pytherm
=======

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
 
  * UNIFAC and Modified UNIFAC 
  * PITZER
  * SIT

* EOS models:

  * PSRK 

* Calculations:

  * Solubility in 2 and 3 component systems 
  * LLE for non-electrolyte systems 
  * VLE using EOS 

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

Roadmap
-------
Future features:

* Activity models:

  * UNIQUAC 
  * NRTL
  * Regular solution
  * LIFAC
  * LIQUAC

* Calculations:

  * Multicomponent equilibrium calculations using PITZER ans SIT models
