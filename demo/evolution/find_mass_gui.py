#!/usr/bin/env python
import sys
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
from matplotlib.widgets import Slider
import stella.evolution as evolution

def find_mass():

    Mbol_sun = 4.74 # according to the XXIX IAU Resolution B2

    find = False
    infile = open('params.dat')
    title = [s.strip() for s in infile.readline().strip()[1:].split(',')]
    for row in infile:
        if row[0] in '%#':
            continue
        g = row.split('|')
        name = g[0].strip().lower().replace(' ','')
        if sys.argv[1].lower()!=name:
            continue
        find = True
        Teff     = int(g[title.index('Teff')])
        Vmag     = float(g[title.index('Vmag')])
        FeH      = float(g[title.index('FeH')])
        Parallax = float(g[title.index('Parallax')])
        BC       = float(g[title.index('BC')])
        break
    infile.close()
    if not find:
        print 'cannot find star %s'%sys.argv[1]
        exit()
    mass    = 1.0
    param   = [mass,FeH]

    Mbol = Vmag + BC - 15. + 5.*math.log10(Parallax) + 5.
    logL = 0.4*(Mbol_sun-Mbol)
    logTeff = math.log10(Teff)

    fig = plt.figure(figsize=(8,6))
    ax  = fig.add_axes([0.15,0.2,0.75,0.70])
    axs = fig.add_axes([0.70,0.25,0.15,0.15])
    ax1 = fig.add_axes([0.15,0.09,0.75,0.03])
    ax2 = fig.add_axes([0.15,0.05,0.75,0.03])
    mbar = Slider(ax1,'$M_*/M_\odot$',0.5,2.5,valinit=mass)
    fbar = Slider(ax2,'$\mathrm{[Fe/H]}$',FeH-1,FeH+1,valinit=FeH)

    def update_mass(val):
        param[0] = val
        replot()
    def update_feh(val):
        param[1] = val
        replot()
    def replot():
        ax.cla()
        track = evolution.get_track('Y2',mass=param[0],feh=param[1])
        logg = math.log10(param[0]) + 4.*math.log10(Teff/5777.)+0.4*(Mbol-4.75)+4.438
        ax.plot(track.logTeff_lst,track.logL_lst,'b-')
        ax.plot(logTeff,logL,'ro',alpha=0.8)
        ax.set_xlim(4.1+1e-4,3.4-1e-4)
        ax.set_ylim(-1.0,4.0)
        ax.set_xlabel('$\log{T_\mathrm{eff}}$')
        ax.set_ylabel('$\log_{10}(L/L_\odot)$')
        ax.text(4.05,3.60,'$M=%4.2f\,M_\odot$'%param[0])
        ax.text(4.05,3.30,'$\mathrm{[Fe/H]}=%+4.2f$'%param[1])
        ax.text(4.05,3.00,'$\log{g}=%4.2f$'%logg)
        ax.xaxis.set_major_locator(tck.MultipleLocator(0.1))
        ax.xaxis.set_minor_locator(tck.MultipleLocator(0.01))
        ax.yaxis.set_major_locator(tck.MultipleLocator(1))
        ax.yaxis.set_minor_locator(tck.MultipleLocator(0.1))

        axs.cla()
        axs.plot(track.logTeff_lst,track.logL_lst,'b-')
        axs.plot(logTeff,logL,'ro',alpha=0.8)
        axs.set_xlim(logTeff+0.01,logTeff-0.01)
        axs.set_ylim(logL-0.05,logL+0.05)
        axs.xaxis.set_major_formatter(tck.FormatStrFormatter('%g'))
        axs.xaxis.set_major_locator(tck.MultipleLocator(0.01))
        axs.xaxis.set_minor_locator(tck.MultipleLocator(0.001))
        axs.yaxis.set_major_locator(tck.MultipleLocator(0.05))
        axs.yaxis.set_minor_locator(tck.MultipleLocator(0.01))

    mbar.on_changed(update_mass)
    fbar.on_changed(update_feh)

    replot()

    plt.show()
    plt.close()

if __name__=='__main__':
    find_mass()
