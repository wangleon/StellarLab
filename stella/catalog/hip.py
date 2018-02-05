import os
import math
from ..utils.fitsio import get_bintable_info
from ..utils.asciitable import structitem_to_dict
from .base import _get_HIP_number

class _HIP(object):
    '''Class for Hipparcos Catalogue (`I/239
    <http://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=I/239>`_, Perryman+
    1997).

    The HIP catalogue contains 118,218 records, with HIP numbers ranging
    from 1 to 120,416.
    For more details, see :ref:`Hipparcos Catalogue<catalog_hip>`.
    '''

    def __init__(self):
        self.catfile = os.path.join(os.getenv('STELLA_DATA'), 'catalog/HIP.fits')
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

    def find_object(self, name, epoch=2000.0, output='dict'):
        '''
        Find record for an object in Hipparcos Catalogue.
    
        .. csv-table:: Descriptions of returned parameters
            :header: Key, Type, Unit, Description
            :widths: 30, 30, 30, 100
    
            HIP,      integer32, ,       HIP number
            RAdeg,    float64,   deg,    Right ascension (*α*) in ICRS at given epoch
            DEdeg,    float64,   deg,    Declination (*δ*) in ICRS at at given epoch
            Vmag,     float32,   mag,    *V* magnitude in Johnson system
            Plx,      float32,   mas,    Parallax
            e_Plx,    float32,   mas,    Error in parallax
            pmRA,     float32,   mas/yr, Proper motion in Right ascension with cos(*δ*) factor
            pmDE,     float32,   mas/yr, Proper motion in Declination
            e_pmRA,   float32,   mas/yr, Error in proper motion in RA
            e_pmDE,   float32,   mas/yr, Error in proper motion in Dec
            BTmag,    float32,   mag,    Mean *B*:sub:`T` magnitude
            e_BTmag,  float32,   mag,    Error in *B*:sub:`T` magnitude
            VTmag,    float32,   mag,    Mean *V*:sub:`T` magnitude
            e_VTmag,  float32,   mag,    Error in *V*:sub:`T` magnitude
            B-V,      float32,   mag,    *B* − *V* color in Johnson system
            e_B-V,    float32,   mag,    Error in *B* − *V* color
            r_B-V,    character, ,       Source of *B* − *V* color
            V-I,      float32,   mag,    *V* − *I* color in Cousin system
            e_V-I,    float32,   mag,    Error in *V* − *I* color
            r_V-I,    character, ,       Source of *V* − *I* color
            Hpmag,    float32,   mag,    Median magnitude in Hipparcos system
            e_Hpmag,  float32,   mag,    Error in *Hp* magnitude
            Hpscat,   float32,   mag,    Scatter on *Hp* magnitude
            o_Hpmag,  integer16, ,       Number of observations for *Hp* magnitude
            SpType,   string,    mag,    Spectral type
            r_SpType, character, ,       Source of Spectral type
    
        Args:
            name (string or integer): Name or number of star.
            epoch (float): Epoch of output astrometric parameters.
            output (string): Type of output results. Either *"dict"* or
                *"dtype"* (:class:`numpy.dtype`).
        Returns:
            dict or :class:`numpy.dtype`: Record in catalogue.
        Raises: 
            FileNotExist: Catalogue file does not exist
        Examples:
            Find record for τ Ceti (HIP 8102)
    
            .. code-block:: python
    
                from stella.catalog import HIP
    
                # find the parametres for tau Cet (HIP 8102)
                item = HIP.find_object(8102, epoch=1991.25)
    
        '''

        def change_epoch(item, epoch):
            pm_ra = item['pmRA']*1e-3/3600. # convert pm_RA from mas/yr to deg/yr
            pm_de = item['pmDE']*1e-3/3600. # convert pm_DE from mas/yr to deg/yr
            item['RAdeg'] += (epoch-1991.25)*pm_ra/math.cos(item['DEdeg']/180.*math.pi)
            item['DEdeg'] += (epoch-1991.25)*pm_de

        hip = _get_HIP_number(name)

        if self._data_info is None:
            self._get_data_info()

        infile = open(self.catfile, 'rb')
        if hip is None:
            # return a null result
            # hip = 672 is the common null record in both HIP and HIP New
            infile.seek(self.pos+(672-1)*self.nbyte, 0)
            item = self.fmtfunc(infile.read(self.nbyte))
        else:
            infile.seek(self.pos+(hip-1)*self.nbyte, 0)
            item = self.fmtfunc(infile.read(self.nbyte))
            change_epoch(item, epoch)
        infile.close()

        if output == 'dtype':
            return item
        elif output == 'dict':
            return structitem_to_dict(item)
        else:
            return None


HIP = _HIP()
