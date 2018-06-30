#!/usr/bin/env python3
import os
import math
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from stella.parameter.teff import _BV_to_Teff_Flower1996
from stella.parameter.bc import _Teff_to_BC_Flower1996
from stella.catalog.utils import plot_skymap

def main():

    # read data file
    catfile = os.path.join(os.getenv('STELLA_DATA'), 'catalog/HIP.fits')
    data = fits.getdata(catfile)
    mask = data['HIP']>0
    data = data[mask]

    # plot skymap
    ra   = data['RAdeg']
    dec  = data['DEdeg']
    plot_skymap(ra, dec, 'skymap_hip.png', size=1, alpha=0.2)

    # plot magnitude histogram
    mask = np.isnan(data['Vmag'])
    vmag = data[~mask]['Vmag']
    bins = np.arange(-2,16)
    fig = plt.figure(figsize=(8,6), dpi=150)
    ax = fig.add_axes([0.1,0.1,0.88,0.85])
    ax.hist(vmag, bins=bins, color='#1166aa', rwidth=0.9)
    ax.set_axisbelow(True)
    ax.set_facecolor('#dddddd')
    ax.yaxis.grid(True, color='w', linestyle='-', linewidth=1)
    ax.set_yscale('log')
    # change tick size
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(13)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(13)
    ax.set_xlabel('$V$', fontsize=15)
    ax.set_ylabel('$N$', fontsize=15)
    ax.set_xticklabels(np.arange(-2, 16, 2))
    # save the figure
    fig.savefig('maghist_hip.png')
    plt.close(fig)

    #---------------------------------------------------------------------------
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

    #---------------------------------------------------------------------------
    fig = plt.figure(figsize=(8,6), dpi=150)
    ax1 = fig.add_axes([0.1,0.1,0.75,0.85])
    ax2 = fig.add_axes([0.88,0.1,0.03,0.85])

    x1, x2, dx = 3000, 12000, 100
    y1, y2, dy = -2, 4, 0.1
    xbins = np.arange(x1, x2+1e-6, dx)
    ybins = np.arange(y1, y2+1e-6, dy)
    nx = xbins - 1
    ny = ybins - 1
    
    norm = mcolors.LogNorm()

    _,_,_,cax = ax1.hist2d(Teff, logL, bins=(xbins, ybins), cmap='Blues', norm=norm)
    cbar = fig.colorbar(cax, cax=ax2)

    # reverse x axis
    _x1, _x2 = ax1.get_xlim()
    ax1.set_xlim(_x2, _x1)

    ax1.set_xlabel('$T_\mathrm{eff}$ (K)', fontsize=15)
    ax1.set_ylabel('$\log(L/L_\odot)$', fontsize=15)
    cbar.set_label('$N$', fontsize=15)

    c1, c2 = cbar.get_clim()
    mticks = []
    power1 = math.floor(math.log10(c1))
    power2 = math.ceil(math.log10(c2))
    for power in np.arange(power1, power2):
        for v in np.arange(1,10)*10**power:
            if c1 <= v <= c2:
                mticks.append(v)
    mticks = cax.norm(mticks)
    cbar.ax.yaxis.set_ticks(mticks, minor=True)

    # change tick size
    for tick in ax1.xaxis.get_major_ticks():
        tick.label1.set_fontsize(13)
    for tick in ax1.yaxis.get_major_ticks():
        tick.label1.set_fontsize(13)
    cbar.ax.tick_params(labelsize=13)

    fig.savefig('hrdhist_hip.png')
    plt.close(fig)


if __name__=='__main__':
    main()
