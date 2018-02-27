import os
import numpy as np
import astropy.io.fits as fits
from ..utils.fitsio import get_bintable_info
from ..utils.asciitable import structitem_to_dict
from .name import _get_HR_number

class _BSC(object):
    '''Class for *Bright Star Catalogue* 5th Edition (`V/50
    <http://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=V/50>`_, Hoffleit+
    1991).

    .. csv-table:: Descriptions of Columns in Catalogue
        :header: Key, Type, Unit, Description
        :widths: 30, 30, 30, 120

        HR,        integer32, ,       HR number
        RAdeg,     float64,   deg,    Right ascension (*α*) in equinox B1950 at epoch 1950.0
        DEdeg,     float64,   deg,    Declination (*δ*) in equinox B1950 at epoch 1950.0
        pmRA,      float32,   mas/yr, Proper motion in Right ascension
        pmDE,      float32,   mas/yr, Proper motion in Declination
        Plx,       float32,   mas,    Parallax
        n_Plx,     character, mas,    Flag of Parallax type
        Vmag,      float32,   mag,    *V* magnitude
        n_Vmag,    character, ,       Code for *V* magnitude
        u_Vmag,    character, ,       Uncertainty flag on *V* magnitude
        B-V,       float32,   mag,    *B* − *V* color in *UBV* system
        u_B-V,     character, ,       Uncertainty flag on *B* − *V*
        U-B,       float32,   mag,    *U* − *B* color in *UBV* system
        u_U-B,     character, ,       Uncertainty flag on *U* − *B*
        R-I,       float32,   mag,    *R* − *I* color
        n_R-I,     character, ,       Code for *R* − *I* color (Cousin or Eggen)
        SpType,    string20,  ,       Spectral type
        n_SpType,  character, ,       Code for spectral type
        RadVel,    float32,   km/s,   Heliocentric radial velocity
        n_RadVel,  string4,   ,       comment on radial velocity
        l_RotVel,  character, ,       limit character on rotational velocity
        RotVel,    float32,   km/s,   Rotational velocity
        u_RotVel,  character, ,       Uncertainty flag on rotation velocity
        Dmag,      float32,   mag,    Magnitude difference of multiple stars
        Sep,       float32,   arcsec, Seperation of components in binary
        MultID,    string4,   ,       Identifications of components in Dmag
        MultCnt,   integer16, ,       Number of components assigned to a multiple
        
    '''

    def __init__(self):
        self.catfile = os.path.join(os.getenv('STELLA_DATA'), 'catalog/BSC.fits')
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
        Find record for an object in *Bright Star Catalogue*, 5th Edition.

        Args:
            name (string or integer): Name or number of star.
            output (string): Type of output results. Either *"dict"* or
                *"dtype"* (:class:`numpy.dtype`).
        Returns:
            dict or :class:`numpy.dtype`: Record in catalogue.
        Examples:

        '''

        hr = _get_HR_number(name)

        if self._data_info is None:
            self._get_data_info()

        nrow    = self._data_info['nrow']
        nbyte   = self._data_info['nbyte']
        pos     = self._data_info['pos']
        fmtfunc = self._data_info['fmtfunc']

        infile = open(self.catfile, 'rb')

        if sao > 0 and hr <= nrow:
            infile.seek(pos+(hr-1)*nbyte,0)
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

BSC = _BSC()
