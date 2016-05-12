#!/usr/bin/env python
import sys,math
import numpy as np
import matplotlib.pyplot as plt
import stella.evolution as evolution

def main():
    trackname = sys.argv[1]
    mass_lst = [1.0,2.0,3.0,4.0]
    feh_lst = [0.0,-0.2,-0.5,-0.8]
    colors = 'rgbc'
    fig = plt.figure(figsize=(12,8))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)
    for imass,mass in enumerate(mass_lst):
        for ifeh,feh in enumerate(feh_lst):
            track = evolution.get_track(trackname,mass=mass,feh=feh)
            if imass == 0:
                ax1.plot(track.logTeff_lst,track.logL_lst,'-',
                    color=colors[ifeh],label='[Fe/H]=%+3.2f'%feh)
                ax2.plot(track.age_lst,track.logL_lst,'-',
                    color=colors[ifeh],label='[Fe/H]=%+3.2f'%feh)
                ax3.plot(track.age_lst,track.logTeff_lst,'-',
                    color=colors[ifeh],label='[Fe/H]=%+3.2f'%feh)
                ax4.plot(track.age_lst,track.R_lst,'-',
                    color=colors[ifeh],label='[Fe/H]=%+3.2f'%feh)
            else:
                ax1.plot(track.logTeff_lst,track.logL_lst,'-',
                    color=colors[ifeh])
                ax2.plot(track.age_lst,track.logL_lst,'-',
                    color=colors[ifeh])
                ax3.plot(track.age_lst,track.logTeff_lst,'-',
                    color=colors[ifeh])
                ax4.plot(track.age_lst,track.R_lst,'-',
                    color=colors[ifeh])
    t1,t2 = ax1.get_xlim()
    l1,l2 = ax1.get_ylim()
    logteffs = np.arange(t1,t2+1e-4,0.1)
    for R in [1,10,100]:
        logl = (math.log10(R) + 2*(logteffs - math.log10(5777.)))*2
        ax1.plot(logteffs,logl,'--',color='gray')
    ax1.set_xlim(t2,t1)
    ax1.set_ylim(l1,l2)


    ax1.set_xlabel('$\log{T_\mathrm{eff}}\,(\mathrm{K})$')
    ax1.set_ylabel('$\log{\left(L/L_\odot\\right)}$')
    ax2.set_xlabel('age (Gyr)')
    ax2.set_ylabel('$\log{\left(L/L_\odot\\right)}$')
    ax3.set_xlabel('age (Gyr)')
    ax3.set_ylabel('$\log{T_\mathrm{eff}}\,(\mathrm{K})$')
    ax4.set_xlabel('age (Gyr)')
    ax4.set_ylabel('$R/R_\odot$')
    leg = ax1.legend(loc='lower left')
    for t in leg.get_texts():
        t.set_fontsize(10)
    leg.get_frame().set_alpha(0.1)
    plt.show()


if __name__=='__main__':
    main()
