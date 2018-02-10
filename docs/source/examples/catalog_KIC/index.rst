.. |Teff| replace:: *T*:sub:`eff`

Kepler Input Catalog
====================

Convert ASCII Files to FITS Table
---------------------------------
Below example converts the ASCII files of KIC to FITS file.

.. literalinclude:: convert.py
   :language: python

Plotting Histograms
-------------------

Below example shows the histogram of *K*:sub:`p` magnitudes and 2D histogram
of KIC.

.. literalinclude:: plot_histogram.py
   :language: python

.. figure:: histogram.png
    :alt: KIC statistics
    :align: center
    :width: 900px
    :figwidth: 900px
    
    Histogram of *K*:sub:`p` magnitudes (*left*) and |Teff| - log\ *g*
    diagram (*right*)

Counting Parameters
-------------------
Below example read the ASCII files of KIC, and counts the number of parallax,
proper motion, and atmospheric parameters in KIC.

.. literalinclude:: count_parameters.py
   :language: python

Results::

    N(PM)=12944973, N(Plx)=958, N(para)=2106821
