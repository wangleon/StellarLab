
.. |EBV| replace:: *E*:sub:`B−V`
.. |AV| replace:: *A*:sub:`V`

Extinction
==========

Dust map
--------
Schlegel et al. 1998 [#Schlegel1998]_ presents a high-resolution Galactic dust
map (hereafter, SFD map) based on Diffuse Infrared Background Experiment on
borad *Cosmic Backgrouund Explorer* (*COBE*/DIRBE) data and *Infrared
Astronomical Satellite* (*IRAS*) observations.

see also `Galactic DUST Reddening & Extinction
<https://irsa.ipac.caltech.edu/applications/DUST/>`_

To get the extinction value (|EBV|) at the Galactic coordinate (*l*, *b*),

.. code-block:: python
    
    from stella.extinction import SFDMap
    EBV = SFDMap.get_EBV(l,b)

.. currentmodule:: stella.extinction
.. autosummary::
   SFDMap

Some studies suggest that this map seems to systematically overesimate the
reddening.
For example, Arce & Goodman 1999 [#Arce1999]_ shows a factor of 1.3-1.5 in
regions of smooth extinction with |Av| > 0.5 mag.
Bonifacio et al. 2000 [#Bonifacio2000]_ uses a revision of |EBV| in SFD map for
|EBV| > 0.10.
Whereas Meléndez et al. 2006 [#Melendez2006]_ adopts both a −10% linear
correction and a zero-point correction of 0.01 mag.
Schlafly & Finkbeiner 2011 [#Schlafly2011]_ suggests a correction of −14% of
SFD map.

In summary,

* Bonifacio et al. (2000) correction:

  .. code-block:: python

      if EBV > 0.1:
          EBV = 0.10 + 0.65*(EBV - 0.10)

* Meléndez et al. (2006) correction:

  .. code-block:: python

      EBV = 0.9*EBV - 0.01

* Schlafly & Finkbeiner (2011) correction:

  .. code-block:: python

      EBV = 0.86*EBV


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
.. [#Melendez2006] Meléndez et al., 2006, *ApJ*, 642, 1082 :ads:`2006ApJ...642.1082M`
.. [#Olsen1988] Olsen, 1988, *A&A*, 189, 173 :ads:`1988A%26A...189..173O`
.. [#Schlafly2011] Schlafly & Finkbeiner, 2011, *ApJ*, 737, 103 :ads:`2011ApJ...737..103S`
.. [#Schlegel1998] Schlegel et al., 1998, *ApJ*, 500, 525 :ads:`1998ApJ...500..525S`
