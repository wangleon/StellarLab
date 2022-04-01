import os
import numpy as np
import astropy.io.fits as fits

from ..utils.download import get_file

def get_phoenix_hires(teff, logg, z, alpha=0.0):
    datapath = 'thirdpartydata/phoenix/v2.0/HiResFITS/'

    # get wavelength from a single fits file
    filepath = os.path.join(datapath, 'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')
    filename = get_file(filepath)
    wave = fits.getdata(filename)

    # get fluxes
    if z == 0.0:
        zcode = '-0.0'
    else:
        zcode = '{:+4.1f}'.format(z)

    if alpha == 0.0:
        acode = ''
    else:
        acode = '.Alpha={:5.2f}'.format(alpha)

    folder1 = 'PHOENIX-ACES-AGSS-COND-2011/'
    folder = 'Z{}{}'.format(zcode, acode)
    fname = 'lte{:05d}-{:4.2f}{}{}.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'.format(
            teff, logg, zcode, acode)
    filepath = os.path.join(datapath, folder1, folder, fname)
    filename = get_file(filepath)
    flux = fits.getdata(filename)
    return wave, flux


def get_phoenix_r10kspec(teff, logg, z, alpha=0.0):
    datapath = 'thirdpartydata/phoenix/v1.0/MedResFITS/R10000FITS/'

    if z == 0.0:
        zcode = '-0.0'
    else:
        zcode = '{:+4.1f}'.format(z)

    if alpha == 0.0:
        acode = ''
    else:
        acode = '.Alpha={:5.2f}'.format(alpha)

    folder = 'PHOENIX-ACES-AGSS-COND-2011_R10000FITS_Z{}{}'.format(zcode, acode)
    fname = 'lte{:05d}-{:4.2f}{}{}.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'.format(
            teff, logg, zcode, acode)
    filepath = os.path.join(datapath, folder, fname)

    filename = get_file(filepath)

    flux, head = fits.getdata(filename, header=True)
    w0 = head['CRVAL1']
    dw = head['CDELT1']
    n = flux.size
    logwave = w0 + np.arange(n)*dw
    wave = np.exp(logwave)
    return wave, flux
