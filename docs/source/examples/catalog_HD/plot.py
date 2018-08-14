#!/usr/bin/env python3
import os
import numpy as np
import astropy.io.fits as fits

from stella.catalog.utils import plot_skymap, plot_histogram

def main():

    # read data file
    path = os.getenv('STELLA_DATA')
    filename = os.path.join(path, 'catalog/HD.fits')
    data = fits.getdata(filename)

    # plot skymap
    ra  = data['RAdeg']
    dec = data['DEdeg']
    plot_skymap(ra, dec, 'skymap_hd.png', size=1, alpha=0.2)

    # plot magnitude histogram
    mask = np.isnan(data['Ptm'])
    ptm = data['Ptm'][~mask]
    plot_histogram(ptm,
            bins    = np.arange(0,15),
            figfile = 'maghist_hd.png',
            xlabel  = '$V$',
            xticks  = np.arange(0, 15, 2),
            )
    
if __name__=='__main__':
    main()
