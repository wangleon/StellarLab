Stellar Parameters
==================

Effective Temperature
---------------------

Photometric Temperature
^^^^^^^^^^^^^^^^^^^^^^^

.. currentmodule:: stella.parameter.teff
.. autosummary::
   color_to_Teff
   _BV_to_Teff_Flower1996
   _get_dwarf_Teff_Alonso1996
   _get_giant_Teff_Alonso1999
   _get_dwarf_Teff_Ramirez2005
   _get_giant_Teff_Ramirez2005
   _get_dwarf_Teff_Masana2006
   _get_dwarf_Teff_GB2009
   _get_giant_Teff_GB2009
   _get_dwarf_Teff_Onehag2009
   _get_giant_Teff_Onehag2009
   _get_dwarf_Teff_Casagrande2010

Bolometric Correction
---------------------
The bolometric correction (BC) is defined as the quantity to be added to the
magnitude in a specific passband in order to obtain its bolometric magnitude.
For example, for *V* band:

.. math::
    \mathrm{BC}_V = m_\mathrm{bol} - m_V = M_\mathrm{bol} - M_V


.. currentmodule:: stella.parameter.bc
.. autosummary::
   get_BC
   _get_dwarf_BC_Alonso1995
   _get_giant_BC_Alonso1999
   _get_dwarf_BC_Masana2006

