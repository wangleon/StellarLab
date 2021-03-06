.. |Msun| replace:: *M*:sub:`⊙`

Stellar Evolution
=================

Geneva Evolutoin Tracks and Ioschrones
--------------------------------------
The Geneva evluion tracks and isochrones (Lejeune & Schaerer 2001
[#Lejeune2001]_) covers mass range of 0.4 ~ 150 |Msun| and metallicity range of
*Z* = 0.0004 ~ 0.1 ([Fe/H] ≈ −1.68 ~ +0.7).

Input parameters
^^^^^^^^^^^^^^^^

* Opacities: OPAL (high temperature, Rogers & Iglesias 1995 [#Rogers1995]_,
  Iglesias & Rogers 1996 [#Iglesias1996]_), Kurucz 1991 [#Kurucz1991]_ or
  Alexander & Ferguson 1994 [#Alexander1994]_ (low temperature);
* Chemical component: The initial (*Y*, *Z*) compostions from linear chemical
  enrichment law *Y* = *Y*:sub:`P` + (Δ\ *Y*/Δ\ *Z*) *Z*, with *Y*:sub:`P` =
  0.24, and Δ\ *Y*/Δ\ *Z* = 3 (for *Z* ≤ 0.02) or 2.5 (for *Z* > 0.02);
* Mass loss rates: Reimers 1975 [#Reimers1975]_ for RGB (red giant brach) and
  EAGB (early asymptotic giant branch); de Jager et al. 1988 [#deJager1988]_ for
  others;
* Core overshooting: *d*/*H*:sub:`P` = 0.2 for stars with *M* ≥ 1.5 |Msun|.
  Tracks with and without overshooting are provided for stars with *M* = 1.25
  |Msun|.

Evolutionary stages
^^^^^^^^^^^^^^^^^^^

* To the end of carbon burning for stars with *M* ≥ 7 |Msun|;
* To the end of early asymptotic gaint branch (EAGB) for stars with 2 |Msun| ≤
  *M* ≤ 5 |Msun|;
* To the stage before the helium flash for stars with 0.8 |Msun| ≤ 1.7 |Msun|.
  Additional tracks include post-helium flashes, i.e. cover the horizontal brach
  (HB) and EAGB for *Z* = 0.020 and 0.001.
* Pre-main sequence (PMS), thermally pulsing AGB, post-AGB, or white dwarfs
  (WDs) are not included in this grid.

Grid parameters
^^^^^^^^^^^^^^^

.. csv-table:: Grid of stellar evolution tracks
   :header: *Z*, *X*, *Y*, *M*:sub:`0` range
   :widths: 10, 10, 10, 40

   0.001, 0.756, 0.243, 0.80 ~ 120.0 |Msun|
   0.004, 0.744, 0.252, 0.80 ~ 120.0 |Msun|
   0.008, 0.728, 0.264, 0.80 ~ 120.0 |Msun|
   0.020, 0.680, 0.300, 0.80 ~ 120.0 |Msun|
   0.040, 0.620, 0.340, 0.80 ~ 120.0 |Msun|
   0.100, 0.420, 0.480, 0.80 ~ 85.0  |Msun|

Module Summary
^^^^^^^^^^^^^^
.. currentmodule:: stella.evolution.geneva
.. autosummary::
    _Geneva.get_track

Yale-Yonsei Evolution Tracks and Ioschrones
-------------------------------------------
The Yale-Yonsei (Y\ :sup:`2`) evolution tracks (Yi et al. 2003 [#Yi2003]_) and
isochrones (Demarque et al. 2004 [#Demarque2004]_)

Input parameters
^^^^^^^^^^^^^^^^
* Opacities: OPAL (high temperature, Rogers & Iglesias 1995 [#Rogers1995]_,
  Iglesias & Rogers 1996 [#Iglesias1996]_), or Alexander & Ferguson 1994
  [#Alexander1994]_ (low temperature);


YaPSI
-----

References
----------
.. [#Alexander1994] Alexander & Ferguson, 1994, *ApJ*, 437, 879 :ads:`1994ApJ...437..879A`
.. [#Demarque2004] Demarque et al., 2004, *ApJS*, 155, 667 :ads:`2004ApJS..155..667D`
.. [#Lejeune2001] Lejeune & Schaerer, 2001, *A&A*, 366, 538 :ads:`2001A%26A...366..538L`
.. [#deJager1988] de Jager et al., 1988, *A&AS*, 72, 259 :ads:`1988A%26AS...72..259D`
.. [#Kurucz1991] Kurucz, 1991, *ASIC*, 341, 441 :ads:`1991ASIC..341..441K`
.. [#Iglesias1996] Iglesias & Rogers, 1996, *ApJ*, 464, 943 :ads:`1996ApJ...464..943I`
.. [#Reimers1975] Reimers, 1975, *MSRSL*, 8, 369 :ads:`1975MSRSL...8..369R`
.. [#Rogers1995] Rogers & Iglesias, 1995, *ASPC*, 78, 31 :ads:`1995ASPC...78...31R`
.. [#Yi2003] Yi et al., 2003, *ApJS*, 144, 259 :ads:`2003ApJS..144..259Y`
