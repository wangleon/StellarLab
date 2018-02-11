#!/usr/bin/env python3
import os
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import matplotlib.ticker as tck

catfile = os.path.join(os.getenv('STELLA_DATA'), 'catalog/EPIC.fits')
data = fits.getdata(catfile)

label_fontsize=15
tick_fontsize=12
fig = plt.figure(figsize=(8,6), dpi=150)
ax1 = fig.add_axes([0.10,0.10,0.75,0.85])
ax3 = fig.add_axes([0.88,0.10,0.02,0.85])

h, yedge, xedge = np.histogram2d(data['logg'], data['Teff'],
                    bins=[np.arange(0,6.01,0.1), np.arange(2000,12001,100)])
cax = ax1.imshow(np.log10(h), cmap='Blues', interpolation='none', aspect='auto')
cbar = fig.colorbar(cax,cax=ax3)
xticks = np.arange(0, 100+1, 20)
ax1.set_xticks(xticks)
ax1.set_xticklabels(np.int32(xedge[xticks]))
yticks = np.arange(0, 60+1, 10)
ax1.set_yticks(yticks)


ax1.set_yticklabels(['%3.1f'%yedge[t] for t in yticks])
ax1.set_xlim(100,0)
ax1.set_ylim(60,0)
ax1.set_xlabel('$T_\mathrm{eff}$ (K)',fontsize=label_fontsize)
ax1.set_ylabel('$\log{g}$',fontsize=label_fontsize)
cbar.set_label('$\log{N}$',fontsize=label_fontsize)
ax1.xaxis.set_minor_locator(tck.MultipleLocator(5))
ax1.yaxis.set_minor_locator(tck.MultipleLocator(5))
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(tick_fontsize)
for tick in ax1.yaxis.get_major_ticks():
    tick.label1.set_fontsize(tick_fontsize)
for tick in cbar.ax.get_yaxis().get_major_ticks():
    tick.label2.set_fontsize(tick_fontsize)
fig.savefig('histogram_hrd.png')
plt.show()
plt.close(fig)
