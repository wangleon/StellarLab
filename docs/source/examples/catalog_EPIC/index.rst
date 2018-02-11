K2 Ecliptic Plane Input catalog (EPIC)
======================================

Convert ASCII file into FITS table
----------------------------------
.. literalinclude:: convert.py
   :language: python

Plotting HRD
------------
Below example shows the distribution of EPIC stars in HR diagram.

.. literalinclude:: plot_hrdpy
   :language: python

.. figure:: histogram_hrd.png
    :alt: EPIC HRD
    :align: center
    :width: 500px
    :figwidth: 500px

    2D histogram of EPIC stars on *T*:sub:`eff` − log\ *g* plane.
