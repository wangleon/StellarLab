#!/usr/bin/env python3
import os
import numpy as np
import astropy.io.fits as fits
from stella.catalog.utils import plot_skymap, plot_histogram

def main():

    # read data file
    catfile = os.path.join(os.getenv('STELLA_DATA'), 'catalog/BSC.fits')
    data = fits.getdata(catfile)

    # plot skymap
    ra   = data['RAdeg']
    dec  = data['DEdeg']
    plot_skymap(ra, dec, 'skymap_bsc.png', size=2, alpha=0.7)

    # plot magnitude histogram
    mask = np.isnan(data['Vmag'])
    vmag = data[~mask]['Vmag']
    plot_histogram(vmag,
            bins    = np.arange(-2, 10),
            figfile = 'maghist_bsc.png',
            xlabel  = '$V$',
            xticks  = np.arange(-2, 9, 2),
            )

if __name__=='__main__':
    main()
