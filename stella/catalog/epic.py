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

        EPIC,      integer32, ,            EPIC Indentifier
        Objtype,   string8,   ,            Object type ('STAR' or 'EXTENDED')
        Kepflag,   string3,   ,            "Kepler magnitude flag ('gri', 'BV', 'JHK', or 'J')"
        RAdeg,     float64,   deg,         Right ascension (*α*) at J2000
        DEdeg,     float64,   deg,         Declination (*δ*) at J2000
        pmRA,      float32,   mas/yr,      Proper motion in Right ascension with cos(*δ*) factor
        pmDE,      float32,   mas/yr,      Proper motion in Declination
        e_pmRA,    float32,   mas/yr,      Error in proper motion in RA
        e_pmDE,    float32,   mas/yr,      Error in proper motion in Dec
        Plx,       float32,   mas,         Parallax
        e_Plx,     float32,   mas,         Error in parallax
        Bmag,      float32,   mag,         *B* magnitude in Johnson system
        e_Bmag,    float32,   mag,         Error in *B* magnitude
        Vmag,      float32,   mag,         *V* magnitude in Johnson system
        e_Vmag,    float32,   mag,         Error in *V* magnitude
        umag,      float32,   mag,         *u* magnitude in Sloan system
        e_umag,    float32,   mag,         Error in *u* magnitude
        gmag,      float32,   mag,         *g* magnitude in Sloan system
        e_gmag,    float32,   mag,         Error in *g* magnitude
        rmag,      float32,   mag,         *r* magnitude in Sloan system
        e_rmag,    float32,   mag,         Error in *r* magnitude
        imag,      float32,   mag,         *i* magnitude in Sloan system
        e_imag,    float32,   mag,         Error in *i* magnitude
        zmag,      float32,   mag,         *z* magnitude in Sloan system
        e_zmag,    float32,   mag,         Error in *z* magnitude
        Jmag,      float32,   mag,         *J* magnitude in 2MASS
        e_Jmag,    float32,   mag,         Error in *J* magnitude
        Hmag,      float32,   mag,         *H* magnitude in 2MASS
        e_Hmag,    float32,   mag,         Error in *H* magnitude
        Kmag,      float32,   mag,         *K* magnitude in 2MASS
        e_Kmag,    float32,   mag,         Error in *K* magnitude
        W1mag,     float32,   mag,         *W*:sub:`1` magnitude in WISE
        e_W1mag,   float32,   mag,         Error in *W*:sub:`1` magnitude
        W2mag,     float32,   mag,         *W*:sub:`2` magnitude in WISE
        e_W2mag,   float32,   mag,         Error in *W*:sub:`2` magnitude
        W3mag,     float32,   mag,         *W*:sub:`3` magnitude in WISE
        e_W3mag,   float32,   mag,         Error in *W*:sub:`3` magnitude
        W4mag,     float32,   mag,         *W*:sub:`4` magnitude in WISE
        e_W4mag,   float32,   mag,         Error in *W*:sub:`4` magnitude
        kepmag,    float32,   mag,         Magnitude in Kepler band 
        Teff,      integer16, K,           Effective temperature
        ue_Teff,   integer16, K,           Upper uncertainty on effective temperature
        le_Teff,   integer16, K,           Lower uncertainty on effective temperature
        logg,      float32,   dex,         Surface gravity
        ue_logg,   float32,   dex,         Upper uncertainty on surface gravity
        le_logg,   float32,   dex,         Lower uncertainty on surface gravity
        FeH,       float32,   dex,         Metallicity
        ue_FeH,    float32,   dex,         Upper uncertainty on metallicity
        le_FeH,    float32,   dex,         Lower uncertainty on metallicity
        Rad,       float32,   *R*:sub:`⊙`, Stellar radius
        ue_Rad,    float32,   *R*:sub:`⊙`, Upper uncertainty on stellar radius
        le_Rad,    float32,   *R*:sub:`⊙`, Lower uncertainty on stellar radius
        Mass,      float32,   *M*:sub:`⊙`, Stellar mass
        ue_Mass,   float32,   *M*:sub:`⊙`, Upper uncertainty on stellar mass
        le_Mass,   float32,   *M*:sub:`⊙`, Lower uncertainty on stellar mass
        rho,       float32,   *ρ*:sub:`⊙`, Stellar density
        ue_rho,    float32,   *ρ*:sub:`⊙`, Upper uncertainty on stellar density
        le_rho,    float32,   *ρ*:sub:`⊙`, Lower uncertainty on stellar density
        Lum,       float32,   *ρ*:sub:`⊙`, Stellar luminosity
        ue_Lum,    float32,   *ρ*:sub:`⊙`, Upper uncertainty on stellar luminosity
        le_Lum,    float32,   *ρ*:sub:`⊙`, Lower uncertainty on stellar luminosity
        Dist,      float32,   pc,          Distance
        ue_Dist,   float32,   pc,          Upper uncertainty on distance
        le_Dist,   float32,   pc,          Lower uncertainty on distance
        E(B-V),    float32,   mag,         Reddening in *B* − *V*
        ue_E(B-V), float32,   mag,         Upper uncertainty on *E*\ (*B* − *V*)
        le_E(B-V), float32,   mag,         Lower uncertainty on *E*\ (*B* − *V*)

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
        '''Get information of FITS table.
        
        Args:
            dataset (integer): Number of EPIC table (1~6).
        Returns:
            No returns

        '''
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

            Find stellar parameters of EPIC 201121245.

            .. code-block:: python

                >>> from stella.catalog import EPIC
                >>> EPIC.find_object(201121245)
                >>> rec = EPIC.find_object(201121245)
                >>> rec['kepmag'], rec['Teff'], rec['logg'], rec['FeH']
                (11.112000465393066, 4899, 2.9719998836517334, -0.3970000147819519)

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
