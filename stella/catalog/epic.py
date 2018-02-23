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

    .. csv-table:: Descriptions of Columns in Catalogue
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

    '''

    def __init__(self):
        self.catfile = {dataset: os.path.join(os.getenv('STELLA_DATA'),
                                             'catalog/EPIC_%d.fits'%dataset)
                        for dataset in range(1,7)
                        }
        self._epic_ranges = {
                1: (201000001, 210000000),
                2: (210000001, 220000000),
                3: (220000001, 230000000),
                4: (230000001, 240000000),
                5: (240000001, 250000000),
                6: (250000001, 251809654),
                }
        self._data_info = {}
        
    def _get_data_info(self, dataset):
        '''Get information of FITS table.'''
        nbyte, nrow, ncol, pos, dtype, fmtfunc = get_bintable_info(self.catfile[dataset])
        self._data_info[dataset] = {
                'nbyte'  : nbyte,
                'nrow'   : nrow,
                'ncol'   : ncol,
                'pos'    : pos,
                'dtype'  : dtype,
                'fmtfunc': fmtfunc
                }
        
    def find_object(self, name, output='dict'):
        '''Find records in *K2 Ecliptic Plane Input Catalog*.

        Args:
            name (string or integer): Name or number of star.
            output (string): Type of output results. Either *"dict"* or
                *"dtype"* (:class:`numpy.dtype`).
        Returns:
            dict or :class:`numpy.dtype`: Record in catalogue.
        Examples:

        '''
        
        epic = _get_EPIC_number(name)
        for dataset, (epic1, epic2) in sorted(self._epic_ranges.items()):
            if epic1 <= epic <= epic2:
                break

        if dataset not in self._data_info:
            self._get_data_info(dataset)

        pos     = self._data_info[dataset]['pos']
        nrow    = self._data_info[dataset]['nrow']
        nbyte   = self._data_info[dataset]['nbyte']
        fmtfunc = self._data_info[dataset]['fmtfunc']

        infile = open(self.catfile[dataset], 'rb')
        infile.seek(pos + (epic-epic1)*nbyte, 0)
        item = fmtfunc(infile.read(nbyte))
        infile.close()

        if item['EPIC'] != epic:
            return None

        if output == 'ndarray':
            return item
        elif output == 'dict':
            return structitem_to_dict(item)
        else:
            return None

EPIC = _EPIC()
