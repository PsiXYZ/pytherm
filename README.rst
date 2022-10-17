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

Current features:
-----------------
* Activity models:
 
  1. UNIFAC and Modified UNIFAC 
  2. PITZER 

* EOS models:

  1. PSRK 

* Calculations:

  1. Solubility in 2 and 3 component systems 
  2. LLE for non-electrolyte systems 
  3. VLE using EOS 

Getting Started
---------------
UNIFAC

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
    print(am.get_y(inp=system, temperature=298))

Some more examples of use can be found in /examples

Roadmap
-------
Future features:

* Activity models:

  1. UNIQUAC 
  2. NRTL 
  3. eNRTL 
  4. LIFAC 

* Calculations:

  1. Multicomponent eqilibrium calcaulations using PITZER 
