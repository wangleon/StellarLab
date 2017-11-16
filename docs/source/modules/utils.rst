Utilities
=========

Summaries
-----------

**Asciifile module (stella.utils.asciifile)**

.. currentmodule:: stella.utils.asciifile

.. autosummary::
    find_sortedfile
    quickfind_sortedfile

**Asciitable module (stella.utils.asciitable)**

.. currentmodule:: stella.utils.asciitable

.. autosummary::
    load_txt
    save_txt

**FITS module (stella.utils.fitsio)**

.. currentmodule:: stella.utils.fitsio

.. autosummary::
    get_bintable_info
    tform_to_dtype
    tform_to_format

**Interpolation module (stella.utils.interpolation)**

.. currentmodule:: stella.utils.interpolation

.. autosummary::
    newton
    parabolic

Functions in utils.asciifile
--------------------------------

.. automodule:: stella.utils.asciifile
    :members:

Functions in utils.asciitable
---------------------------------

.. automodule:: stella.utils.asciitable
    :members:


Functions in utils.fitsio
--------------------------------

.. automodule:: stella.utils.fitsio
    :members:


Functions in utils.interpolation
----------------------------------

.. automodule:: stella.utils.interpolation
    :members:

Visualization
--------------

Temperature to RGB color
^^^^^^^^^^^^^^^^^^^^^^^^^


.. currentmodule:: stella.utils.vision

.. autosummary::
    temp_to_rgb

The above function converts black-body temperature to RGB color using piecewise
polynomials.
The coefficients of the polynomials are obtained by fitting the CIE 1964 10-deg
color matching functions given in `bbr_color.txt
<http://www.vendian.org/mncharity/dir3/blackbody/UnstableURLs/bbr_color.html>`_
on `Mitchell Charity's website <http://www.vendian.org/mncharity/dir3/blackbody/>`_.


.. plot:: plot/temp_to_rgb.py
