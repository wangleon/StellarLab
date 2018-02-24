#!/usr/bin/env python3
import os
import numpy as np
import astropy.io.fits as fits
from stella.catalog.utils import plot_skymap, plot_histogram, plot_histogram2d

def main():

    ra_lst  = np.array([])
    dec_lst = np.array([])
    kp_lst  = np.array([])
    teff_lst = np.array([])
    logg_lst = np.array([])
    for i in range(1,7):
        # read data file
        catfile = os.path.join(os.getenv('STELLA_DATA'), 'catalog/EPIC_%d.fits'%i)
        data = fits.getdata(catfile)
        mask = data['EPIC']>0

        # plot skymap
        ra_lst  = np.append(ra_lst, data[mask]['RAdeg'])
        dec_lst = np.append(dec_lst, data[mask]['DEdeg'])
        kp_lst  = np.append(kp_lst, data[mask]['kepmag'])
        mask1 = data['Teff']>0
        teff_lst = np.append(teff_lst, data[mask1]['Teff'])
        logg_lst = np.append(logg_lst, data[mask1]['logg'])

    plot_skymap(ra_lst, dec_lst, 'skymap_epic.png', size=1, alpha=0.01)

    # plot magnitude histogram
    bins = np.arange(0,20)
    plot_histogram(kp_lst, '$K_\mathrm{p}$', 'maghist_epic.png', bins=bins, yscale='log')

    # plot Kiel histogram
    plot_histogram2d(teff_lst, logg_lst, (2000, 12000, 100), (0,6,0.1),
            xlabel='$T_\mathrm{eff}$ (K)',
            ylabel='$\log{g}$',
            reverse_x=True, reverse_y=True,
            figfile='kielhist_epic.png')

if __name__=='__main__':
    main()
