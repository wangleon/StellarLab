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

References
----------
.. [#Arce1999] Arce & Goodman, 1999, *ApJ*, 512, L135 :ads:`1999ApJ...512L.135A`
.. [#Bonifacio2000] Bonifacio et al., 2000, *AJ*, 120, 2065 :ads:`2000AJ....120.2065B`
.. [#Schlegel1998] Schlegel et al., 1998, *ApJ*, 500, 525 :ads:`1998ApJ...500..525S`
