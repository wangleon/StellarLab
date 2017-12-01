#!/usr/bin/env python3
import os
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as tck
import matplotlib.cm as cmap
from stella.evolution.geneva import read_track

def main():

    feh1, feh2 = -3, 0.5
    cnorm = colors.Normalize(vmin=feh1, vmax=feh2)
    scalarmap = cmap.ScalarMappable(norm=cnorm, cmap=cmap.jet)
    
    family = 'serif'
    fig = plt.figure(figsize=(8,6), dpi=150)
    ax  = fig.add_axes([0.1,0.1,0.7,0.85])
    axc = fig.add_axes([0.84,0.1,0.03,0.85])
    for iz, z in enumerate([0.020, 0.004, 0.001]):
        M_H = math.log10(z/0.0134)
        color = scalarmap.to_rgba(M_H)
        for im,m0 in enumerate([1.0, 2.0, 5.0]):
            track = read_track(mass0=m0, z=z)
            logTeff_lst = track[0]
            logL_lst    = track[1]
            if im == 0:
                ax.plot(logTeff_lst, logL_lst, '-', color=color, zorder=3-iz,
                        label='$Z = %5.3f\ (\mathrm{[Fe/H]} = %+5.2f)$'%(z, M_H))
            else:
                ax.plot(logTeff_lst, logL_lst, '-', color=color, zorder=3-iz)
            if iz == 0:
                ax.text(logTeff_lst[0]+0.1, logL_lst[0]-0.3,
                        '$M = %3.1f M_\odot$'%m0)

    cax = ax.scatter([0], [0], c=[0], norm=cnorm, cmap=scalarmap.get_cmap())
    leg = ax.legend(loc='lower left')
    for t in leg.get_texts():
        t.set_fontsize(10)
    ax.set_xlabel('$\log(T_\mathrm{eff})$')
    ax.set_ylabel('$\log(L/L_\odot)$')
    cbar = fig.colorbar(cax, cax=axc, ticks=np.arange(feh1, feh2+1e-5, 0.5))
    cbar.set_label('[Fe/H]', family=family)
    ax.set_xlim(4.4, 3.5)
    ax.set_ylim(-1, 5)
    ax.xaxis.set_major_locator(tck.MultipleLocator(0.2))
    ax.xaxis.set_minor_locator(tck.MultipleLocator(0.02))
    ax.yaxis.set_major_locator(tck.MultipleLocator(1))
    ax.yaxis.set_minor_locator(tck.MultipleLocator(0.1))
    cbar.set_clim(feh1, feh2)

    ax.grid(True)
    fig.savefig('HRD_Geneva_track_lowmass.png')
    plt.show()

if __name__=='__main__':
    main()
