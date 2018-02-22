Catalog
=======

Base functions
--------------
.. currentmodule:: stella.catalog.base

.. autosummary::
    _str_to_float
    _str_to_int

.. automodule:: stella.catalog.base
    :members:
    :private-members:
    :undoc-members:

Parsing Starnames
-----------------
.. currentmodule:: stella.catalog.name
.. autosummary::
   _get_regular_BD_name
   _get_regular_CD_name
   _get_regular_G_name
   _get_regular_HD_name
   _get_regular_HIP_name
   _get_regular_TYC_name
   _get_star_number1
   _get_EPIC_number
   _get_HIP_number
   _get_KIC_number
   _get_TYC_number


.. automodule:: stella.catalog.name
    :members:
    :private-members:
    :undoc-members:

Catalogues
----------

.. currentmodule:: stella.catalog

.. autosummary::
   hip._HIP
   hip._HIP2
   tyc2._TYC2
   kic._KIC
   epic._EPIC

.. automodule:: stella.catalog.hip
   :members:
   :private-members:
   :undoc-members:

.. automodule:: stella.catalog.tyc2
   :members:
   :private-members:
   :undoc-members:

.. automodule:: stella.catalog.kic
   :members:
   :private-members:
   :undoc-members:

.. automodule:: stella.catalog.epic
   :members:
   :private-members:
   :undoc-members:

.. automodule:: stella.catalog.kepler
   :members:
   :private-members:
   :undoc-members:

Cross Index
------------
.. currentmodule:: stella.catalog.xindex

.. autosummary::
    BD_to_HIP
    CD_to_HIP
    G_to_TYC
    HD_to_HIP
    HD_to_TYC
    HIP_to_2MASS
    HIP_to_Gaia
    HIP_to_HD
    HIP_to_BD
    HIP_to_CD
    HIP_to_TYC
    Kepler_to_KIC
    Kepler_to_KOI
    KIC_to_Kepler
    KIC_to_KOI
    KOI_to_Kepler
    KOI_to_KIC
    TYC_to_2MASS
    TYC_to_HIP

.. automodule:: stella.catalog.xindex
    :members:
    :private-members:
    :undoc-members:

Catalogue Utilities
-------------------

.. currentmodule:: stella.catalog.utils

.. autosummary::
    plot_skymap
    plot_histogram
    plot_histogram2d

.. automodule:: stella.catalog.utils
    :members:
    :private-members:
    :undoc-members:
