from __future__ import print_function
import os
import math
import numpy as np
import astropy.io.fits as fits
from scipy.interpolate import RectBivariateSpline

class SFDMapClass(object):
    """Galactic dust map of `Schlegel+ 1998
    <http://adsabs.harvard.edu/abs/1998ApJ...500..525S>`_.

    Examples:
    .. code-block:: python

        from stella.extinction import SFDMap
        EBV = SFDMap.get_EBV(l,b)

    where *l*, *b* are Galactic coordinate.
    """
    def __init__(self):

        self.data = {'n': None, 's': None}
        self.head = {'n': None, 's': None}

    def _read_data(self, key):
        """Read FITS data.

        Args:
            key (str): Either 'n' (northern hemisphere) or 's' (southern
                hemisphere)
        """
        filename = os.path.join(os.getenv('STELLA_DATA'),
                    'extinction/SFD_dust_4096_%sgp.fits'%key)
        self.data[key], self.head[key] = fits.getdata(filename, header=True)

    def get_EBV(self, l, b):
        """get *E(B-V)* from the SFD dust map.

        Args:
            l (float): Galactic longitude in degree
            b (float): Galactic latitude in degree
        Returns:
            float: E(B-V)

        .. code-block:: python
            
            from stella.extinction import SFDMap
            EBV = SFDMap.get_EBV(l,b)

        """
        if b >= 0:
            key = 'n'
        else:
            key = 's'

        if self.data[key] is None:
            self._read_data(key)

        data = self.data[key]
        head = self.head[key]
        naxis1  = head['NAXIS1']
        naxis2  = head['NAXIS2']
        ctype1  = head['CTYPE1']
        ctype2  = head['CTYPE2']
        crval1  = head['CRVAL1']
        crval2  = head['CRVAL2']
        crpix1  = head['CRPIX1']
        crpix2  = head['CRPIX2']
        lonpole = head['LONPOLE']
        if 'CDELT1' in head and 'CDELT2' in head:
            cd1_1 = head['CDELT1']
            cd2_2 = head['CDELT2']
            cd1_2 = 0
            cd2_1 = 0
        else:
            cd1_1 = head['CD1_1']
            cd1_2 = head['CD1_2']
            cd2_1 = head['CD2_1']
            cd2_2 = head['CD2_2']
    
        if ctype1=='LAMBERT--X' and ctype2=='LAMBERT--Y':
            nsgp  = head['LAM_NSGP']
            scale = head['LAM_SCAL']
            rho = math.sqrt(2)*math.sin((45.-0.5*nsgp*b)/180.*math.pi)
            xpix =         scale * rho * math.cos(l/180.*math.pi)
            ypix = -nsgp * scale * rho * math.sin(l/180.*math.pi)
            xpix = xpix + (crpix1 - crval1 - 1.0)
            ypix = ypix + (crpix2 - crval2 - 1.0)
        elif ctype1[-3:]=='ZEA' and ctype2[-3:]=='ZEA':
            if abs(crval2 - 90) < 1e-3:
                theta = b
                phi = (l + 180. + lonpole - crval1) % 360
            elif abs(crval2 + 90) < 1e-3:
                theta = -b
                phi = (lonpole + crval1 - l) % 360
            rtheta = 2.* 180./math.pi * math.sin(0.5*(90.-theta)/180.*math.pi)
            xtemp =  rtheta * math.sin(phi/180.*math.pi)
            ytemp = -rtheta * math.cos(phi/180.*math.pi)
            if cd1_1 == 0 and cd2_2 == 0:
                xpix = xtemp/cd1_1 + (crpix1 - 1.0)
                ypix = ytemp/cd2_2 + (crpix2 - 1.0)
            else:
                denom = cd1_1 * cd2_2 - cd1_2 * cd2_1
                xpix = (cd2_2 * xtemp - cd1_2 * ytemp) / denom + (crpix1 - 1.)
                ypix = (cd1_1 * ytemp - cd2_1 * xtemp) / denom + (crpix2 - 1.)
        else:
            print('Unsupported projection method in CTYPE keywords')
    
        x1 = min(max(int(xpix)-3,0),naxis1-6)
        x2 = x1 + 6
        y1 = min(max(int(ypix)-3,0),naxis2-6)
        y2 = y1 + 6
        sdata = data[y1:y2,x1:x2]
    
        y = np.arange(data.shape[0])[y1:y2]
        x = np.arange(data.shape[1])[x1:x2]
        fnew = RectBivariateSpline(y, x, sdata)
        EBV = fnew(ypix, xpix)[0][0]

        return EBV

SFDMap = SFDMapClass()
