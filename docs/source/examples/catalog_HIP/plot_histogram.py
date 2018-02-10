#!/usr/bin/env python3
import os
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
from stella.parameter.teff import _BV_to_Teff_Flower1996
from stella.parameter.bc import _Teff_to_BC_Flower1996

catfile = os.path.join(os.getenv('STELLA_DATA'), 'catalog/HIP.fits')
data = fits.getdata(catfile)

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

label_fontsize=15
tick_fontsize=12
bins = np.arange(-2,15)
fig = plt.figure(figsize=(15,6), dpi=150)
ax1 = fig.add_axes([0.07,0.12,0.39,0.83])
ax2 = fig.add_axes([0.52,0.12,0.39,0.83])
ax3 = fig.add_axes([0.92,0.12,0.015,0.83])
ax1.hist(data['Vmag'], bins=bins, color='#1166aa', lw=0.5)
ax1.set_yscale('log')
ax1.set_xlabel('$V$ Magnitude', fontsize=label_fontsize)
ax1.set_ylabel('$N$', fontsize=label_fontsize)

h, yedge, xedge = np.histogram2d(logL, Teff,
                    bins=[np.arange(-2,4+1e-3,0.1), np.arange(3000,12001,100)])
cax = ax2.imshow(np.log10(h), cmap='Blues', interpolation='none', aspect='auto')
cbar = fig.colorbar(cax,cax=ax3)

xticks = np.arange(0, 90+1, 20)
ax2.set_xticks(xticks)
ax2.set_xticklabels(np.int32(xedge[xticks]))
yticks = np.arange(0, 60+1, 10)
ax2.set_yticks(yticks)
#for t in yticks:
#    print(t, yedge[t])
ax2.set_yticklabels(['%3.1f'%yedge[t] for t in yticks])
ax2.set_xlim(90,0)
ax2.set_ylim(0,60)
ax2.set_xlabel('$T_\mathrm{eff}$ (K)',fontsize=label_fontsize)
ax2.set_ylabel('$\log(L/L_\odot)$',fontsize=label_fontsize)
cbar.set_label('$\log{N}$',fontsize=label_fontsize)
#ax2.xaxis.set_major_locator(tck.MultipleLocator(1000))
ax2.xaxis.set_minor_locator(tck.MultipleLocator(5))
#ax2.yaxis.set_major_locator(tck.MultipleLocator(1.0))
ax2.yaxis.set_minor_locator(tck.MultipleLocator(5))

for ax in [ax1, ax2]:
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(tick_fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(tick_fontsize)
for tick in cbar.ax.get_yaxis().get_major_ticks():
    tick.label2.set_fontsize(tick_fontsize)

fig.savefig('histogram.png')
plt.show()
plt.close(fig)
