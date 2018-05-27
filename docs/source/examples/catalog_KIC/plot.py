#!/usr/bin/env python3
import os
import numpy as np
import astropy.io.fits as fits
from stella.catalog.utils import plot_skymap, plot_histogram, plot_histogram2d

def main():

    # read data file
    catfile = os.path.join(os.getenv('STELLA_DATA'), 'catalog/KIC.fits')
    data = fits.getdata(catfile)
    #mask = data['HIP']>0
    #data = data[mask]

    # plot skymap
    ra   = data['RAdeg']
    dec  = data['DEdeg']
    plot_skymap(ra, dec, 'skymap_kic.png', size=1, alpha=0.2)

    exit()

    # plot magnitude histogram
    mask = np.isnan(data['Vmag'])
    vmag = data[~mask]['Vmag']
    bins = np.arange(-2,16)
    plot_histogram(vmag, '$V$', 'maghist_hip.png', bins=bins, yscale='log')

    # plot HRD
    mask1 = ~np.isnan(data['Plx'])
    mask2 = ~np.isnan(data['B-V'])
    data = data[mask1*mask2]
    mask = data['Plx']>0
    data = data[mask]
    Mv = data['Vmag'] - 5*np.log10(1000./data['Plx']) + 5
    Teff = _BV_to_Teff_Flower1996(data['B-V'])
    mask = Teff>0
    Teff=Teff[mask]
    BC = np.array([_Teff_to_BC_Flower1996(t) for t in Teff])
    Mbol = Mv[mask] + BC
    logL = 0.4*(4.74-Mbol)

    plot_histogram2d(Teff, logL, (3000, 12000, 100), (-2,4,0.1),
            xlabel='$T_\mathrm{eff}$ (K)',
            ylabel='$\log(L/L_\odot)$',
            reverse_x=True,
            figfile='hrdhist_hip.png')

if __name__=='__main__':
    main()
