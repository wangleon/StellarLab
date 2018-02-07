Extinction
==========

Dust map
--------
Schlegel et al. 1998 [#Schlegel1998]_ presents Galactic dust map based on COBE
and IRAS data.

also read http://www.astro.princeton.edu/~Schlegel/dust/local/local.html

Use SFD map:

.. code-block:: python
    
    from stella.extinction import SFDMap
    EBV = SFDMap.get_EBV(l,b)

reduce correction by Bonifacio et al. 2000 [#Bonifacio2000]_.
Also read Arce & Goodman 1999 [#Arce1999]_

.. code-block:: python

    if EBV > 0.1:
        EBV = 0.10 + 0.65*(EBV - 0.10)

For stars in reddening layer, correct the E(B-V) by a factor of
:math:`1-\exp(-|D\sin b|/h)`, where *D* is the distance in pc, *b* is the
galactic latitude, and *h* = 125 pc is the scale height of the Galactic
reddening layer.

.. code-block:: python

    import math
    EBV *= 1-math.exp(-abs(d*math.sin(b/180.*math.pi))/125.)

Photometric Extinction
----------------------
Olsen 1988 [#Olsen1988]_ gave the calibration relations of *E*\ (*b* − *y*) v.s.
Strömgren *uvbyβ* colors.
This relation is based on 1231 stars, and valid for stars with spectral types of
F0-G2 and luminosity classses of III-V.

.. currentmodule:: stella.extinction.photometric
.. autosummary::
   get_Stromgren_Eby


References
----------
.. [#Arce1999] Arce & Goodman, 1999, *ApJ*, 512, L135 :ads:`1999ApJ...512L.135A`
.. [#Bonifacio2000] Bonifacio et al., 2000, *AJ*, 120, 2065 :ads:`2000AJ....120.2065B`
.. [#Olsen1988] Olsen, 1988, *A&A*, 189, 173 :ads:`1988A%26A...189..173O`
.. [#Schlegel1998] Schlegel et al., 1998, *ApJ*, 500, 525 :ads:`1998ApJ...500..525S`
