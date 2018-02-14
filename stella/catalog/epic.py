import os
import struct
import numpy as np
import astropy.io.fits as fits
from ..utils.fitsio import get_bintable_info
from ..utils.asciitable import structitem_to_dict
from .name import _get_EPIC_number

class _EPIC(object):
    '''Class for *K2 Ecliptic Plane Input Catalog* (EPIC, `Huber+ 2016
    <http://adsabs.harvard.edu/abs/2016ApJS..224....2H>`_).

    For more details, see :ref:`K2 Ecliptic Plane Input Catalog<catalog_epic>`.
    '''

    def __init__(self):
        self.catfile = os.path.join(os.getenv('STELLA_DATA'), 'catalog/EPIC.fits')
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
        '''Find records in K2 Ecliptic Plane Input Catalog.

        .. csv-table:: Descriptions of returned parameters
            :header: Key, Type, Unit, Description
            :widths: 30, 30, 30, 120

            EPIC,     integer32, ,            EPIC Indentifier
            Teff,     integer16, K,           Effective temperature
            E_Teff,   integer16, K,           Upper uncertainty on effective temperature
            e_Teff,   integer16, K,           Lower uncertainty on effective temperature
            logg,     float32,   dex,         Surface gravity
            E_logg,   float32,   dex,         Upper uncertainty on surface gravity
            e_logg,   float32,   dex,         Lower uncertainty on surface gravity
            FeH,      float32,   dex,         Metallicity
            E_FeH,    float32,   dex,         Upper uncertainty on metallicity
            e_FeH,    float32,   dex,         Lower uncertainty on metallicity
            Rad,      float32,   *R*:sub:`⊙`, Stellar radius
            E_Rad,    float32,   *R*:sub:`⊙`, Upper uncertainty on stellar radius
            e_Rad,    float32,   *R*:sub:`⊙`, Lower uncertainty on stellar radius
            Mass,     float32,   *M*:sub:`⊙`, Stellar mass
            E_Mass,   float32,   *M*:sub:`⊙`, Upper uncertainty on stellar mass
            e_Mass,   float32,   *M*:sub:`⊙`, Lower uncertainty on stellar mass
            rho,      float32,   *ρ*:sub:`⊙`, Stellar density
            E_rho,    float32,   *ρ*:sub:`⊙`, Upper uncertainty on stellar density
            e_rho,    float32,   *ρ*:sub:`⊙`, Lower uncertainty on stellar density
            Dist,     float32,   pc,          Distance
            E_Dist,   float32,   pc,          Upper uncertainty on distance
            e_Dist,   float32,   pc,          Lower uncertainty on distance
            E(B-V),   float32,   mag,         Reddening in *B* − *V*
            E_E(B-V), float32,   mag,         Upper uncertainty on *E*\ (*B* − *V*)
            e_E(B-V), float32,   mag,         Lower uncertainty on *E*\ (*B* − *V*)
            Flag,     string3,   ,            Classification Flag
            RAdeg,    float64,   deg,         Right ascension (*α*) at J2000
            DEdeg,    float64,   deg,         Declination (*δ*) at J2000

        Args:
            name (string or integer): Name or number of star.
            output (string): Type of output results. Either *"dict"* or
                *"dtype"* (:class:`numpy.dtype`).
        Returns:
            dict or :class:`numpy.dtype`: Record in catalogue.
        Examples:

        '''
        
        epic = _get_EPIC_number(name)

        if self._data_info is None:
            self._get_data_info()

        pos     = self._data_info['pos']
        nrow    = self._data_info['nrow']
        nbyte   = self._data_info['nbyte']
        fmtfunc = self._data_info['fmtfunc']

        infile = open(self.catfile, 'rb')

        # the EPIC numbers in the FITS file are increased monotonoically.
        i1, i2 = 0, nrow-1
        while(i2-i1 > 1):
            i3 = int((i1+i2)/2)
            infile.seek(pos + i3*nbyte, 0)
            key = struct.unpack('>i',infile.read(4))[0]
            if epic < key:
                i2 = i3
            elif epic > key:
                i1 = i3
            else:
                break
        infile.seek(-4,1)
        item = fmtfunc(infile.read(nbyte))

        infile.close()

        if output == 'ndarray':
            return item
        elif output == 'dict':
            return structitem_to_dict(item)
        else:
            return None

EPIC = _EPIC()
