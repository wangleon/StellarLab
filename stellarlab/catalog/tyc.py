import os
import math
import struct
import numpy as np
from ..utils.fitsio import get_bintable_info
from ..utils.asciitable import structitem_to_dict
from .name import _get_TYC_number

class _TYC(object):
    """Class for *Tycho Catalogue* (`I/239/tyc_main
    <http://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=I/239/tyc_main>`_, ESA
    1997).

    The *Tycho Catalogue* contains 1,058,332 records.

    .. csv-table:: Descriptions of Columns in Catalogue
        :header: Key, Type, Unit, Description
        :widths: 30, 30, 30, 180

        TYC,      integer32, ,       TYC number
        RAdeg,    float64,   deg,    Right ascension (*α*) in ICRS at given epoch
        DEdeg,    float64,   deg,    Declination (*δ*) in ICRS at at given epoch
        e_RA,     float32,   mas,    Error in Right ascension
        e_DE,     float32,   mas,    Error in Declination
        pmRA,     float32,   mas/yr, Proper motion in Right ascension with cos(*δ*) factor
        pmDE,     float32,   mas/yr, Proper motion in Declination
        e_pmRA,   float32,   mas/yr, Error in proper motion in RA
        e_pmDE,   float32,   mas/yr, Error in proper motion in Dec
        Plx,      float32,   mas,    Parallax
        e_Plx,    float32,   mas,    Error in parallax
        Nastro,   integer16, ,       Number of transits for astrometry
        Vmag,     float32,   mag,    *V* magnitude in Johnson system
        r_Vmag,   character, ,       Source of magnitude
        BTmag,    float32,   mag,    Mean *B*:sub:`T` magnitude
        e_BTmag,  float32,   mag,    Error in *B*:sub:`T` magnitude
        VTmag,    float32,   mag,    Mean *V*:sub:`T` magnitude
        e_VTmag,  float32,   mag,    Error in *V*:sub:`T` magnitude
        r_BTmag,  character, ,       Source of photometry
        B-V,      float32,   mag,    *B* − *V* color in Johnson system
        e_B-V,    float32,   mag,    Error in *B* − *V* color
        Nphoto,   integer16, ,       Number of transits for photometry
        VTscat,   float32,   mag,    Scatter of *V*:sub:`T` magnitude
        VTmax,    float32,   mag,    maximum *V*:sub:`T` magnitude (15%)
        VTmin,    float32,   mag,    minimum *V*:sub:`T` magnitude (85%)
        VarFlag,  character, ,       Variability flag
        MultFlag, character, ,       Duplicity flag
        
    """

    def __init__(self):
        self.catfile = os.path.join(os.getenv('STELLA_DATA'), 'catalog/TYC.fits')
        self._data_info = None

    def _get_data_info(self):
        """Get information of FITS table."""
        nbyte, nrow, ncol, pos, dtype, fmtfunc = get_bintable_info(self.catfile)
        self._data_info = {
                'nbyte'  : nbyte,
                'nrow'   : nrow,
                'ncol'   : ncol,
                'pos'    : pos,
                'dtype'  : dtype,
                'fmtfunc': fmtfunc
                }

    def find_object(self, name, epoch=2000.0, output='dict'):
        """Find record for an object in *Tycho Catalogue*.
    
        Args:
            name (string or tuple): Name or nunmber tuple of star.
            epoch (float): Epoch of output astrometric parameters.
            output (string): Type of output results. Either *"dict"* or
                *"dtype"* (:class:`numpy.dtype`).
        Returns:
            dict or :class:`numpy.dtype`: Record in catalogue. If input name
                can't be find in this catalogue, an empty item with TYC=0 is
                returned.
        Examples:
            Find proper motion of AD Leo (TYC 1423-174-1)
    
            .. code-block:: python
            
                >>> from stella.catalog import TYC
                >>> rec = TYC.find_object('TYC 1423-174-1')
                >>> rec['pmRA'], rec['pmDE']
                (-516.4000244140625, -25.299999237060547)

        """

        tyc1, tyc2, tyc3 = _get_TYC_number(name)
    
        target = np.int32((tyc1<<18) + (tyc2<<4) + (tyc3<<1))

        if self._data_info is None:
            self._get_data_info()

        pos     = self._data_info['pos']
        nbyte   = self._data_info['nbyte']
        nrow    = self._data_info['nrow']
        fmtfunc = self._data_info['fmtfunc']
    
        infile = open(self.catfile, 'rb')
        i1, i2 = 1, nrow-2
        find = False
        while(i2-i1 > 1):
            i3 = (i1+i2)//2
            infile.seek(pos + i3*nbyte, 0)
            key = struct.unpack('>i',infile.read(4))[0]
            if target < key:
                i2 = i3
            elif target > key:
                i1 = i3
            else:
                find = True
                break

        if find:
            infile.seek(-4,1)
        else:
            # the first element is an empty item
            infile.seek(pos, 0)
        item = fmtfunc(infile.read(nbyte))

        infile.close()

        # change epoch
        pm_ra = item['pmRA']*1e-3/3600. # convert pm_RA from mas/yr to deg/yr
        pm_de = item['pmDE']*1e-3/3600. # convert pm_DE from mas/yr to deg/yr
        item['RAdeg'] += (epoch-1991.25)*pm_ra/math.cos(item['DEdeg']/180.*math.pi)
        item['DEdeg'] += (epoch-1991.25)*pm_de
    
        if output == 'ndarray':
            return item
        elif output == 'dict':
            return structitem_to_dict(item)
        else:
            return None

class _TYC2(object):
    """Class for *Tycho-2 Catalogue* (`I/259
    <http://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=I/259>`_, Høg+ 2000).

    The *Tycho-2 Catalogue* contains 2,539,913 records, and additonal 17,588
    records in Supplement-1, and 1,146 records in Supplement-2.
    For more details, see :ref:`Tycho-2 Catalogue<catalog_tyc2>`.

    .. csv-table:: Descriptions of Columns in Catalogue
        :header: Key, Type, Unit, Description
        :widths: 30, 30, 30, 180
    
        TYC,      integer32, ,       TYC number
        RAdeg,    float64,   deg,    Right ascension (*α*) in ICRS at given epoch
        DEdeg,    float64,   deg,    Declination (*δ*) in ICRS at at given epoch
        e_RA,     integer16, mas,    Error in RA with cos(*δ*) factor at mean epoch (−1 if blank)
        e_DE,     integer16, mas,    Error in Dec at mean epoch (−1 if blank)
        epoch_RA, float32,   yr,     Mean epoch of RA in Julian year
        epoch_DE, float32,   yr,     Mean epoch of Dec in Julian year
        q_RA,     float32,   ,       Goodness of fit for mean RA (truncated to 9.9 if >9.9)
        q_DE,     float32,   ,       Goodness of fit for mean Dec (truncated to 9.9 if >9.9)
        pmRA,     float32,   mas/yr, Proper motion in Right ascension with cos(*δ*) factor
        pmDE,     float32,   mas/yr, Proper motion in Declination
        e_pmRA,   float32,   mas/yr, Error in proper motion in RA
        e_pmDE,   float32,   mas/yr, Error in proper motion in Dec
        q_pmRA,   float32,   ,       Goodness of fit for pmRA (truncated to 9.9 if >9.9)
        q_pmDE,   float32,   ,       Goodness of fit for pmDE (truncated to 9.9 if >9.9)
        Num,      integer16, ,       Number of positions used (−1 if blank)
        BTmag,    float32,   mag,    Mean *B*:sub:`T` magnitude
        e_BTmag,  float32,   mag,    Error in *B*:sub:`T` magnitude
        VTmag,    float32,   mag,    Mean *V*:sub:`T` magnitude
        e_VTmag,  float32,   mag,    Error in *V*:sub:`T` magnitude
        prox,     float32,   arcsec, "Proximity, or distance to the nearest entry (truncated to 99.9 if >99.9)"

    """

    def __init__(self):
        self.catfile = os.path.join(os.getenv('STELLA_DATA'), 'catalog/TYC2.fits')
        self._data_info = None

    def _get_data_info(self):
        """Get information of FITS table."""
        nbyte, nrow, ncol, pos, dtype, fmtfunc = get_bintable_info(self.catfile)
        self._data_info = {
                'nbyte'  : nbyte,
                'nrow'   : nrow,
                'ncol'   : ncol,
                'pos'    : pos,
                'dtype'  : dtype,
                'fmtfunc': fmtfunc
                }

    def find_object(self, name, epoch=2000.0, output='dict'):
        """Find record for an object in *Tycho-2 Catalogue*.
    
        Args:
            name (string or tuple): Name or nunmber tuple of star.
            epoch (float): Epoch of output astrometric parameters.
            output (string): Type of output results. Either *"dict"* or
                *"dtype"* (:class:`numpy.dtype`).
        Returns:
            dict or :class:`numpy.dtype`: Record in catalogue. If input name
                can't be find in this catalogue, an empty item with TYC=0 is
                returned.
        Examples:
            Find proper motion of Barnard's star (TYC 425-2502-1)
    
            .. code-block:: python
            
                >>> from stella.catalog import TYC2
                >>> rec = TYC2.find_object('TYC 425-2502-1')
                >>> rec['pmRA'], rec['pmDE']
                (-798.7999877929688, 10277.2998046875)

        """

        tyc1, tyc2, tyc3 = _get_TYC_number(name)
    
        target = np.int32((tyc1<<18) + (tyc2<<4) + (tyc3<<1))

        if self._data_info is None:
            self._get_data_info()

        pos     = self._data_info['pos']
        nbyte   = self._data_info['nbyte']
        nrow    = self._data_info['nrow']
        fmtfunc = self._data_info['fmtfunc']
    
        infile = open(self.catfile, 'rb')
        i1, i2 = 1, nrow-2
        find = False
        while(i2-i1 > 1):
            i3 = (i1+i2)//2
            infile.seek(pos + i3*nbyte, 0)
            key = struct.unpack('>i',infile.read(4))[0]
            if target < key:
                i2 = i3
            elif target > key:
                i1 = i3
            else:
                find = True
                break

        if find:
            infile.seek(-4,1)
        else:
            # the first element is an empty item
            infile.seek(pos, 0)
        item = fmtfunc(infile.read(nbyte))
        
        # change epoch
        pm_ra = item['pmRA']*1e-3/3600. # convert pm_RA from mas/yr to deg/yr
        pm_de = item['pmDE']*1e-3/3600. # convert pm_DE from mas/yr to deg/yr
        item['RAdeg'] += (epoch-2000.0)*pm_ra/math.cos(item['DEdeg']/180.*math.pi)
        item['DEdeg'] += (epoch-2000.0)*pm_de

        # looking for possible companion
        if i3 != nrow - 1:
            key2 = struct.unpack('>i',infile.read(4))[0]
            if key2 == key + 1:
                infile.seek(-4,1)
                item2 = fmtfunc(infile.read(nbyte))
                print('Warning: There are more than 1 star matched')
                t1 = (key2 & 0b11111111111111000000000000000000)/2**18
                t2 = (key2 & 0b00000000000000111111111111110000)/2**4
                t3 = (key2 & 0b00000000000000000000000000001110)/2
                print(t1, t2, t3, item2)
    
        infile.close()
    
        if output == 'ndarray':
            return item
        elif output == 'dict':
            return structitem_to_dict(item)
        else:
            return None

TYC2 = _TYC2()
TYC  = _TYC()
