#!/usr/bin/env python3
import os
import numpy as np
import astropy.io.fits as fits
from stella.catalog.utils import plot_skymap, plot_histogram

def main():

    # read data file
    catfile = os.path.join(os.getenv('STELLA_DATA'), 'catalog/SAO.fits')
    data = fits.getdata(catfile)

    # plot skymap
    ra   = data['RAdeg']
    dec  = data['DEdeg']
    plot_skymap(ra, dec, 'skymap_sao.png', size=1, alpha=0.2)

    # plot magnitude histogram
    mask = np.isnan(data['Vmag'])
    vmag = data[~mask]['Vmag']
    bins = np.arange(0,14)
    plot_histogram(vmag, '$V$', 'maghist_sao.png', bins=bins, yscale='log')


if __name__=='__main__':
    main()
