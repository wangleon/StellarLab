import os
import numpy as np
import astropy.io.fits as fits
from ..utils.fitsio import get_bintable_info
from ..utils.asciitable import structitem_to_dict
from .name import _get_HD_number

class _HD(object):
    '''Class for *Hennry Draper Catalogue* (`III/135A
    <http://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=III/135A>`_, Cannon &
    Pickering 1918-1924).

    .. csv-table:: Descriptions of Columns in Catalogue
        :header: Key, Type, Unit, Description
        :widths: 30, 30, 30, 120

        HD,    integer32, ,    HD number
        RAdeg, float64,   deg, Right ascension (*α*) in equinox B1900 at epoch 1900.0
        DEdeg, float64,   deg, Declination (*δ*) in equinox B1900 at epoch 1900.0
        q_Ptm, integer16, ,    "Code for Ptm: 0 = measured, 1 = inferred from Ptg and spectral type"
        Ptm,   float32,   mag, Photovisual magnitude
        q_Ptg, integer16, ,    "Code for Ptg: 0 = measured, 1 = inferred from Ptm and spectral type"
        Ptg,   float32,   mag, Photographic magnitude
        SpT,   string,    ,    Spectral type
        Int,   string,    ,    Photographic intensity of spectrum
        Rem,   character, ,    Remarks

    '''

    def __init__(self):
        self.catfile = os.path.join(os.getenv('STELLA_DATA'), 'catalog/HD.fits')
        self._data_info = None

    def _get_data_info(self):
        '''Get information of FITS table.'''
        nbyte, nrow, ncol, pos, dtype, fmtfunc = get_bintable_info(self.catfile)
        self._data_info = {
                'nbyte'  : nbyte,
                'nrow'   : nrow,
                'ncol'   : ncol,
                'pos'    : pos,
                'dtype'  : dtype,
                'fmtfunc': fmtfunc
                }

    def find_object(self, name, output='dict'):
        '''
        Find record for an object in *Henry Draper Catalogue*.

        Args:
            name (string or integer): Name or number of star.
            output (string): Type of output results. Either *"dict"* or
                *"dtype"* (:class:`numpy.dtype`).
        Returns:
            dict or :class:`numpy.dtype`: Record in catalogue.
        Examples:
            Find record for τ Ceti (HD 10700)

            .. code-block:: python

                >>> from stella.catalog import HD
                >>> rec = HD.find_object('HD 10700')
                >>> rec['RAdeg'], rec['DEdeg'], rec['Ptm'], rec['Ptg'], rec['SpT']
                (24.85, -16.466666666666665, 3.6500000953674316, 4.650000095367432, 'K0')

        '''

        hd = _get_HD_number(name)

        if self._data_info is None:
            self._get_data_info()

        nrow    = self._data_info['nrow']
        nbyte   = self._data_info['nbyte']
        pos     = self._data_info['pos']
        fmtfunc = self._data_info['fmtfunc']

        infile = open(self.catfile, 'rb')

        if hd > 0 and hd <= nrow:
            infile.seek(pos+(hd-1)*nbyte,0)
            item = fmtfunc(infile.read(nbyte))
        else:
            pass

        infile.close()

        if output == 'ndarray':
            return item
        elif output == 'dict':
            return structitem_to_dict(item)
        else:
            return None

HD = _HD()
