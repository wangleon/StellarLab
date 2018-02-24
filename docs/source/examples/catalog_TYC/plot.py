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
    bins = np.arange(0,16)
    plot_histogram(vtmag, '$V_\mathrm{T}$', 'maghist_tyc2.png', bins=bins, yscale='log')


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
    bins = np.arange(0,16)
    plot_histogram(vtmag, '$V_\mathrm{T}$', 'maghist_tyc.png', bins=bins, yscale='log')


    exit()
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
