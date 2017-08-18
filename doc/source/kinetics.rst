
.. |Msun| replace:: *M*:sub:`⊙`
.. |kms| replace:: km s\ :sup:`−1`

Kinetics
=========

Summaries
---------

**Stellar kinetics**

.. currentmodule:: stella.kinetics.orbit
.. autosummary::
    compute_UVW
    compute_GalXYZ


**Galactic potentials**

.. currentmodule:: stella.kinetics.potential
.. autosummary::
    Potential
    PointPotential
    HernquistPotential
    MiyamotoNagaiPotential
    NFWPotential

Stellar orbits in the Milky Way
--------------------------------
Stellar orbits in the Milky Way can be calculated by numerially integrating the
equations of motion.
The Galactic potenial can be approximated with the sum of three components
(bulge, disk, and halo, e.g., Gnedin et al. 2005 [#Gnedin2005]_),

.. math::

   \Phi = \Phi_\mathrm{bluge} + \Phi_\mathrm{disk} + \Phi_\mathrm{halo}

or four components (including the central black hole, e.g., Kenyon et al. 2008
[#Kenyon2008]_).

.. math::

   \Phi = \Phi_\mathrm{BH} + \Phi_\mathrm{bluge} + \Phi_\mathrm{disk} + \Phi_\mathrm{halo}

Disk Potential
^^^^^^^^^^^^^^^

The disk potential can be described by the asymmetric Miyamoto & Nagai 1975
[#Miyamoto1975]_ model:

.. math::

    \Phi_\mathrm{disk}(R, z) = -\frac{GM_\mathrm{disk}}{\sqrt{R^2 + (a+\sqrt{z^2 + b^2})^2}}

where (*R*, *z*) is the radial distance and the height in the galactocentric
cylindrical coordinate system, *M*:sub:`disk` is mass of disk, *a* is the scale
length and *b* is the scale height. The table below lists different parameter
adopted in various studies.

======================== ========== ========== ==============================
*M*:sub:`disk` (|Msun|)  *a* (kpc)  *b* (kpc)  References
======================== ========== ========== ==============================
5.0 × 10\ :sup:`10`      2.4         0.18      Dehnen & Binney 1998 [#Dehnen1998a]_
4.0 × 10\ :sup:`10`      5           0.3       Gnedin et al. 2005 [#Gnedin2005]_
4.0 × 10\ :sup:`10`      5           0.3       Yu & Madau 2007 [#Yu2007]_
6 × 10\ :sup:`10`        2.75        0.3       Kenyon et al. 2014 [#Kenyon2014]_
5.36 × 10\ :sup:`10`     2.75        0.3       Fragione et al. 2017 [#Fragione2017]_
======================== ========== ========== ==============================

Orbital Parameters
^^^^^^^^^^^^^^^^^^^

Serveral orbital parameters are defined to characterize stellar orbits in the
Milky Way (e.g. Takeda 2007 [#Takeda2007]_):

* The mean galactocentric radius *R*:sub:`mean` ≡ (*R*:sub:`max` + *R*:sub:`min`)/2;

* The ellipticity *e* ≡ (*R*:sub:`max` − *R*:sub:`min`)/(*R*:sub:`max` + *R*:sub:`min`);

* The maximum separation from the Galactic plane (*z*:sub:`max`);

* The absolute speed relative to the local standard of rest (\|\ *v*\ \|\ :sub:`LSR`);

References
------------
.. [#Dehnen1998a] Dehnen & Binney, 1998, *MNRAS*, 294, 429 :ads:`1998MNRAS.294..429D`
.. [#Fragione2017] Fragione et al., 2017, *MNRAS*, 467, 451 :ads:`2017MNRAS.467..451F`
.. [#Gnedin2005] Gnedin et al., 2005, *ApJ*, 634, 344 :ads:`2005ApJ...634..344G`
.. [#Kenyon2008] Kenyon et al., 2008, *ApJ*, 680, 312 :ads:`2008ApJ...680..312K`
.. [#Kenyon2014] Kenyon et al., 2014, *ApJ*, 793, 122 :ads:`2014ApJ...793..122K`
.. [#Miyamoto1975] Miyamoto & Nagai, 1975, *PASJ*, 27, 533 :ads:`1975PASJ...27..533M`
.. [#Takeda2007] Takeda, 2007, *PASJ*, 59, 335 :ads:`2007PASJ...59..335T`
.. [#Yu2007] Yu & Madau, 2007, *MNRAS*, 379, 1293 :ads:`2007MNRAS.379.1293Y`
