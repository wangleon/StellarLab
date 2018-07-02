#!/usr/bin/env python3
import os
import math
import numpy as np
import astropy.io.fits as fits
from stella.catalog.utils import plot_skymap, plot_histogram, plot_histogram2d

def main():

    # read data file
    catfile = os.path.join(os.getenv('STELLA_DATA'), 'catalog/KIC.fits')
    data = fits.getdata(catfile)

    # plot skymap
    ra   = data['RAdeg']
    dec  = data['DEdeg']
    plot_skymap(ra, dec, 'skymap_kic.png', size=1, alpha=0.2)

    # plot magnitude histogram
    mask = np.isnan(data['kepmag'])
    kpmag = data['kepmag'][~mask]

    plot_histogram(kpmag,
            bins    = np.arange(0, 26),
            figfile = 'maghist_kic.png',
            xlabel  ='$K_\mathrm{p}$ Magnitude',
            xticks  = np.arange(0, 27, 2),
            )

    #---------------------------------------------------------------------------
    # plot HRD
    mask1 = np.isnan(data['Teff'])
    mask2 = np.isnan(data['logg'])
    mask3 = (~mask1)*(~mask2)
    data3 = data[mask3]
    Teff = data3['Teff']
    logg = data3['logg']

    plot_histogram2d(Teff, logg,
            xbins     = np.arange(3000, 12001, 100),
            ybins     = np.arange(0, 6.01, 0.1),
            xlabel    = '$T_\mathrm{eff}$ (K)',
            ylabel    = '$\log{g}$',
            figfile   = 'kielhist_kic.png',
            reverse_x = True,
            reverse_y = True,
            scale     = 'log',
            )
            
if __name__=='__main__':
    main()
