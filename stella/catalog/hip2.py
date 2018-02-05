import os
import math
from ..utils.fitsio import get_bintable_info
from ..utils.asciitable import structitem_to_dict
from .base import _get_HIP_number

class _HIP2(object):
    '''Class for Hipparcos New Reduction (`I/311
    <http://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=I/311>`_, van Leeuwen
    2007).

    The Hipparcos Catalogue New Reduction contains 117,955 records, with HIP
    numbers ranging from 1 to 120,404.
    For more details, see :ref:`Hipparcos Catalogue New Reduction<catalog_hip2>`.
    '''

    def __init__(self):
        self.catfile = os.path.join(os.getenv('STELLA_DATA'), 'catalog/HIP2.fits')
        nbyte, nrow, ncol, pos, dtype, fmtfunc = get_bintable_info(self.catfile)
        self.nbyte   = nbyte
        self.pos     = pos
        self.fmtfunc = fmtfunc

    def find_object(self, name, epoch=2000.0, output='dict'):
        '''
        Find record for an object in Hipparcos Catalogue New Reduction.

        .. csv-table:: Descriptions of returned parameters
            :header: Key, Type, Unit, Description
            :widths: 30, 30, 30, 100
    
            HIP,     integer32, ,       HIP number
            RAdeg,   float64,   deg,    Right ascension (*α*) in ICRS at given epoch
            DEdeg,   float64,   deg,    Declination (*δ*) in ICRS at at given epoch
            Plx,     float32,   mas,    Parallax
            e_Plx,   float32,   mas,    Error in parallax
            pmRA,    float32,   mas/yr, Proper motion in Right ascension with cos(*δ*) factor
            pmDE,    float32,   mas/yr, Proper motion in Declination
            e_pmRA,  float32,   mas/yr, Error in proper motion in RA
            e_pmDE,  float32,   mas/yr, Error in proper motion in Dec
            B-V,     float32,   mag,    *B* − *V* color in Johnson system
            e_B-V,   float32,   mag,    Error in *B* − *V* color
            V-I,     float32,   mag,    *V* − *I* color in Cousin system
            Hpmag,   float32,   mag,    Median magnitude in Hipparcos system
            e_Hpmag, float32,   mag,    Error in *Hp* magnitude
            Hpscat,  float32,   mag,    Scatter on *Hp* magnitude
    
        Args:
            name (string or integer): Name or number of star.
            epoch (float): Epoch of output astrometric parameters.
            output (string): Type of output results. Either *"dict"* or *"dtype"*
                (:class:`numpy.dtype`).
        Returns:
            dict or :class:`numpy.dtype`: Record in catalogue.
        Raises:
            FileNotExist: Catalogue file does not exist
        Examples:
            Find the ICRS coordinate of τ Ceti (HIP 8102)
    
            .. code-block:: python
            
                from stella.catalog import HIP2
        
                res1 = HIP2.find_object(8102, epoch=1991.25)
                res2 = HIP2.find_object(8102)
                print(res1['RAdeg'], res1['DEdeg'])
                print(res2['RAdeg'], res2['DEdeg'])
                # output:
                # 26.021364586713265 -15.939555724635493
                # 26.017014215650022 -15.937479641367434
        '''

        def change_epoch(item, epoch):
            pm_ra = item['pmRA']*1e-3/3600. # convert pm_RA from mas/yr to deg/yr
            pm_de = item['pmDE']*1e-3/3600. # convert pm_DE from mas/yr to deg/yr
            item['RAdeg'] += (epoch-1991.25)*pm_ra/math.cos(item['DEdeg']/180.*math.pi)
            item['DEdeg'] += (epoch-1991.25)*pm_de

        hip = _get_HIP_number(name)

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

HIP2 = _HIP2()
