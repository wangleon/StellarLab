#!/usr/bin/env python3
import os
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt

from stella.parameter.teff import _BV_to_Teff_Flower1996
from stella.parameter.bc import _Teff_to_BC_Flower1996
from stella.catalog.utils import plot_skymap


def plot_histogram(vtmag, figfile):
    # plot magnitude histogram
    bins = np.arange(-2,16)
    fig = plt.figure(figsize=(8,6), dpi=150)
    ax = fig.add_axes([0.1,0.1,0.88,0.85])
    ax.hist(vtmag, bins=bins, color='#1166aa', rwidth=0.9)
    ax.set_axisbelow(True)
    ax.set_facecolor('#dddddd')
    ax.yaxis.grid(True, color='w', linestyle='-', linewidth=1)
    ax.set_yscale('log')
    # change tick size
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(13)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(13)
    ax.set_xlabel('$V_\mathrm{T}$', fontsize=15)
    ax.set_ylabel('$N$', fontsize=15)
    ax.set_xticklabels(np.arange(-2, 16, 2))
    ax.set_ylim(0.5, 2e6)
    # save the figure
    fig.savefig(figfile)
    plt.close(fig)

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
    plot_histogram(vtmag, 'maghist_tyc2.png')

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
    plot_histogram(vtmag, 'maghist_tyc.png')

if __name__=='__main__':
    main()
