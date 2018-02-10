#!/usr/bin/env python3
import os
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import matplotlib.ticker as tck

catfile = os.path.join(os.getenv('STELLA_DATA'), 'catalog/KIC.fits')

data = fits.getdata(catfile)
mask = np.isnan(data['kepmag'])
kps = data['kepmag'][~mask]

mask1 = np.isnan(data['Teff'])
mask2 = np.isnan(data['logg'])
mask3 = (~mask1)*(~mask2)
data3 = data[mask3]

family = 'Sans Serif'
label_fontsize=15
tick_fontsize=12
bins = np.arange(0,26)
fig1 = plt.figure(figsize=(15,6), dpi=150)
ax1 = fig1.add_axes([0.07,0.12,0.39,0.83])
ax2 = fig1.add_axes([0.52,0.12,0.39,0.83])
ax3 = fig1.add_axes([0.92,0.12,0.015,0.83])
ax1.hist(kps, bins=bins, color='#1166aa', lw=0.5)
ax1.set_yscale('log')
ax1.set_xlabel('$K_\mathrm{p}$ Magnitude', fontsize=label_fontsize,family=family)
ax1.set_ylabel('$N$', fontsize=label_fontsize,family=family)

h, yedge, xedge = np.histogram2d(data3['logg'], data3['Teff'],
                    bins=[np.arange(0,6.01,0.1), np.arange(3000,12001,100)])
cax = ax2.imshow(np.log10(h), cmap='Blues', interpolation='none', aspect='auto')
cbar = fig1.colorbar(cax,cax=ax3)
xticks = np.arange(0, 90+1, 20)
ax2.set_xticks(xticks)
ax2.set_xticklabels(np.int32(xedge[xticks]))
yticks = np.arange(0, 60+1, 10)
ax2.set_yticks(yticks)
#for t in yticks:
#    print(t, yedge[t])
ax2.set_yticklabels(['%3.1f'%yedge[t] for t in yticks])
ax2.set_xlim(90,0)
ax2.set_ylim(60,0)
ax2.set_xlabel('$T_\mathrm{eff}$ (K)',fontsize=label_fontsize,family=family)
ax2.set_ylabel('$\log{g}$',fontsize=label_fontsize,family=family)
cbar.set_label('$\log{N}$',fontsize=label_fontsize,family=family)
#ax2.xaxis.set_major_locator(tck.MultipleLocator(1000))
ax2.xaxis.set_minor_locator(tck.MultipleLocator(5))
#ax2.yaxis.set_major_locator(tck.MultipleLocator(1.0))
ax2.yaxis.set_minor_locator(tck.MultipleLocator(5))

for ax in [ax1, ax2]:
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(tick_fontsize)
        tick.label1.set_family(family)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(tick_fontsize)
        tick.label1.set_family(family)
for tick in cbar.ax.get_yaxis().get_major_ticks():
    tick.label2.set_fontsize(tick_fontsize)
    tick.label2.set_family(family)
fig1.savefig('histogram.png')
plt.show()
plt.close(fig1)
