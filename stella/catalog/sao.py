import os
import numpy as np
import astropy.io.fits as fits
from ..utils.fitsio import get_bintable_info
from ..utils.asciitable import structitem_to_dict
from .name import _get_SAO_number

class _SAO(object):
    '''Class for SAO catalog (`I/131A <>`_).

    .. csv-table:: Descriptions of Columns in Catalogue
        :header: Key, Type, Unit, Description
        :widths: 30, 30, 30, 120

        SAO,       integer32, ,       SAO number
        RAdeg,     float64,   deg,    Right ascension (*α*) in equinox B1950 at epoch 1950.0
        DEdeg,     float64,   deg,    Declination (*δ*) in equinox B1950 at epoch 1950.0
        pmRA,      float32,   mas/yr, Proper motion in Right ascension
        pmDE,      float32,   mas/yr, Proper motion in Declination
        e_pmRA,    float32,   mas/yr, Error in proper motion in RA
        e_pmDE,    float32,   mas/yr, Error in proper motion in Dec
        Vmag,      float32,   mag,    *V* magnitude
        Pmag,      float32,   mag,    Photographic magnitude
        SpType,    string3,   ,       Spectral type
        r_Vmag,    integer16, ,       Source of *V* magnitude
        r_Pmag,    integer16, ,       Source of Photographic magnitude
        r_Num,     integer16, ,       Source of star number
        r_SpType,  integer16, ,       Source of spectral type
        RAdeg2000, float64,   deg,    Right ascension (*α*) in equinox J2000.0 at epoch J2000.0
        DEdeg2000, float64,   deg,    Declination (*δ*) in equinox J2000.0 at epoch J2000.0
        pmRA2000,  float32,   mas/yr, Proper motion in Right ascension in FK5 system
        pmDE2000,  float32,   mas/yr, Proper motion in Declination in FK5 system

    '''

    def __init__(self):
        self.catfile = os.path.join(os.getenv('STELLA_DATA'), 'catalog/SAO.fits')
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
                'fmtfunc': fmtfunc,
                }

    def find_object(self, name, output='dict'):
        '''
        Find record for an object in *SAO Catalogue*.

        Args:
            name (string or integer): Name or number of star.
            output (string): Type of output results. Either *"dict"* or
                *"dtype"* (:class:`numpy.dtype`).
        Returns:
            dict or :class:`numpy.dtype`: Record in catalogue.
        Examples:
            Find record for τ Ceti (SAO 147986)

            .. code-block:: python

                >>> from stella.catalog import SAO
                >>> rec = SAO.find_object('SAO 147986')
                >>> rec['RAdeg2000'], rec['DEdeg2000']
                (26.017158333333334, -15.937394444444445)

        '''

        sao = _get_SAO_number(name)

        if self._data_info is None:
            self._get_data_info()

        nrow    = self._data_info['nrow']
        nbyte   = self._data_info['nbyte']
        pos     = self._data_info['pos']
        fmtfunc = self._data_info['fmtfunc']

        infile = open(self.catfile, 'rb')

        if sao > 0 and sao <= nrow:
            infile.seek(pos+(sao-1)*nbyte,0)
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

SAO = _SAO()
