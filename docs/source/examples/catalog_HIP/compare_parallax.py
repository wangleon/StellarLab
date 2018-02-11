#!/usr/bin/env python3
import os
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt

catfile1 = os.path.join(os.getenv('STELLA_DATA'), 'catalog/HIP.fits')
catfile2 = os.path.join(os.getenv('STELLA_DATA'), 'catalog/HIP2.fits')

data1 = fits.getdata(catfile1)
data2 = fits.getdata(catfile2)
imax = min(data1.size, data2.size)
data1 = data1[0:imax]
data2 = data2[0:imax]

mask1 = ~np.isnan(data1['Plx'])
mask2 = ~np.isnan(data2['Plx'])
mask = mask1*mask2

data1 = data1[mask]
data2 = data2[mask]

mask1 = data1['Plx']>0
mask2 = data2['Plx']>0
mask = mask1*mask2

data1 = data1[mask]
data2 = data2[mask]

m = ((data1['Plx']+data2['Plx'])/2.)<5

fig = plt.figure(figsize=(15,5),dpi=150)
ax1 = fig.add_axes([0.06, 0.13, 0.27, 0.81])
ax2 = fig.add_axes([0.38, 0.13, 0.27, 0.81])
ax3 = fig.add_axes([0.70, 0.13, 0.27, 0.81])
ax1.plot(data1['Plx'][m], data2['Plx'][m], 'ro', alpha=0.1, ms=1)
ax1.plot(data1['Plx'][~m], data2['Plx'][~m], 'go', alpha=0.1, ms=1)
#ax2.plot(data1['e_Plx']/data1['Plx'], data2['e_Plx']/data2['Plx'], 'ko', alpha=0.1, ms=1)
re1 = data1['e_Plx']/data1['Plx']
re2 = data2['e_Plx']/data2['Plx']
ax2.plot(re1[~m], re2[~m], 'go', alpha=0.1, ms=1)
ax3.plot(re1[m], re2[m], 'ro', alpha=0.1, ms=1)
for ax in [ax1, ax2, ax3]:
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True)
    x1, x2 = ax.get_xlim()
    y1, y2 = ax.get_ylim()
    z1 = min(x1, y1)
    z2 = max(x2, y2)
    if ax == ax1:
        z1, z2 = z2, z1
    else:
        z1, z2 = 1e-4, 1e3
    ax.plot([z1, z2], [z1, z2], '-', lw=0.5, color='#1166aa')
    ax.set_xlim(z1, z2)
    ax.set_ylim(z1, z2)
ax1.set_xlabel('$\\varpi$ (mas) - Perryman 1997')
ax1.set_ylabel('$\\varpi$ (mas) - van Leeuwen 2007')
ax2.set_xlabel('$\Delta\\varpi/\\varpi$ (Perryman 1997)')
ax2.set_ylabel('$\Delta\\varpi/\\varpi$ (van Leeuwen 2007)')
ax3.set_xlabel('$\Delta\\varpi/\\varpi$ (Perryman 1997)')
ax3.set_ylabel('$\Delta\\varpi/\\varpi$ (van Leeuwen 2007)')
fig.savefig('compare_parallax.png')
plt.show()
