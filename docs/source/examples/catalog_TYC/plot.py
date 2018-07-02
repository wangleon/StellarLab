#!/usr/bin/env python3
import os
import numpy as np
import astropy.io.fits as fits

from stella.parameter.teff import _BV_to_Teff_Flower1996
from stella.parameter.bc import _Teff_to_BC_Flower1996
from stella.catalog.utils import plot_skymap, plot_histogram, plot_histogram2d

def main():

    # read data file
    catfile = os.path.join(os.getenv('STELLA_DATA'), 'catalog/TYC2.fits')
    data = fits.getdata(catfile)

    # plot skymap
    ra   = data['RAdeg']
    dec  = data['DEdeg']
    plot_skymap(ra, dec, 'skymap_tyc2.png', size=1, alpha=0.1)

    # plot magnitude histogram
    mask = np.isnan(data['VTmag'])
    vtmag = data[~mask]['VTmag']
    plot_histogram(vtmag,
            bins    = np.arange(-2,16),
            figfile = 'maghist_tyc2.png',
            xlabel  = '$V_\mathrm{T}$',
            ylim    = (0.5, 2e6),
            xticks  = np.arange(-2, 17, 2))

    #---------------------------------------------------------------------------
    # read data file
    catfile = os.path.join(os.getenv('STELLA_DATA'), 'catalog/TYC.fits')
    data = fits.getdata(catfile)

    # plot skymap
    ra   = data['RAdeg']
    dec  = data['DEdeg']
    plot_skymap(ra, dec, 'skymap_tyc.png', size=1, alpha=0.1)

    # plot magnitude histogram
    mask = np.isnan(data['VTmag'])
    vtmag = data[~mask]['VTmag']
    plot_histogram(vtmag,
            bins    = np.arange(-2,16),
            figfile = 'maghist_tyc.png',
            xlabel  = '$V_\mathrm{T}$',
            ylim    = (0.5, 2e6),
            xticks  = np.arange(-2, 17, 2))

if __name__=='__main__':
    main()
