#!/usr/bin/env python3
import os
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt

for dataset in ['HIP', 'HIP2']:
    catfile = os.path.join(os.getenv('STELLA_DATA'), 'catalog/%s.fits'%dataset)

    data = fits.getdata(catfile)
    mask = ~np.isnan(data['Plx'])
    data = data[mask]
    mask = data['Plx']>0
    data = data[mask]
    
    fig = plt.figure(figsize=(14,6), dpi=150)
    ax1 = fig.add_axes([0.08, 0.1, 0.27, 0.85])
    ax2 = fig.add_axes([0.38, 0.1, 0.27, 0.85])
    ax3 = fig.add_axes([0.68, 0.1, 0.27, 0.85])
    para = data['Plx']
    d    = 1000./para
    re   = data['e_Plx']/para
    mag  = data['Hpmag']
    ax1.plot(para,re,'ko', ms=1, alpha=0.1)
    ax2.plot(d,re,   'ko', ms=1, alpha=0.1)
    ax3.plot(mag,re, 'ko', ms=1, alpha=0.1)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax3.set_yscale('log')
    for ax in fig.get_axes():
        x1,x2 = ax.get_xlim()
        #y1,y2 = ax.get_ylim()
        #ax.plot([x1,x2],[0.1,0.1],'r-')
        #ax.plot([x1,x2],[0.2,0.2],'r-')
        ax.set_xlim(x1,x2)
        ax.set_ylim(1e-4,1e3)
        ax.grid(True)
    
    ax2.set_yticklabels([])
    ax3.set_yticklabels([])
    
    ax1.set_ylabel('$\Delta\\varpi/\\varpi$')
    ax1.set_xlabel('$\\varpi$ (mas)')
    ax2.set_xlabel('Distance (pc)')
    ax3.set_xlabel('$H_\mathrm{p}$ magnitude')
    fig.savefig('parallax_error_%s.png'%dataset)
    plt.close(fig)

