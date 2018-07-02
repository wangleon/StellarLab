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

    plot_histogram(kp_lst,
            bins    = np.arange(0, 20),
            figfile = 'maghist_epic.png',
            xlabel  = '$K_\mathrm{p}$',
            xticks  = np.arange(0, 21, 2),
            )

    # plot Kiel histogram
    plot_histogram2d(teff_lst, logg_lst,
            xbins     = np.arange(2000, 12001, 100),
            ybins     = np.arange(0, 6.01, 0.1),
            xlabel    = '$T_\mathrm{eff}$ (K)',
            ylabel    = '$\log{g}$',
            figfile   = 'kielhist_epic.png',
            reverse_x = True,
            reverse_y = True,
            scale     = 'log',
            )
    exit()

    # plot magnitude histogram
    bins = np.arange(0,20)
    fig = plt.figure(figsize=(8,6), dpi=150)
    ax = fig.add_axes([0.1,0.1,0.88,0.85])
    ax.hist(kp_lst, bins=bins, color='#1166aa', rwidth=0.9)
    ax.set_axisbelow(True)
    ax.set_facecolor('#dddddd')
    ax.yaxis.grid(True, color='w', linestyle='-', linewidth=1)
    ax.set_yscale('log')
    # change tick size
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(13)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(13)
    ax.set_xlabel('$K_\mathrm{p}$', fontsize=15)
    ax.set_ylabel('$N$', fontsize=15)
    ax.set_xticks(np.arange(0, 20+1e-3, 2))
    ax.set_xlim(0, 20)
    # save the figure
    fig.savefig('maghist_epic.png')
    plt.close(fig)




if __name__=='__main__':
    main()
