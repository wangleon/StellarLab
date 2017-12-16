#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from stella.parameter.teff import color_to_Teff

def main():

    BVtry_lst = np.arange(0.0,1.8+1e-4,0.01)
    ref_lst = ['Alonso1996','Ramirez2005','GB2009','Casagrande2010']

    fig = plt.figure(figsize=(15,5), dpi=150)
    for ifeh, feh in enumerate([-3,-2,-1,0]):
        BV_lst, teff_lst = {}, {}
        for ref in ref_lst:
            BV_lst[ref] = []
            teff_lst[ref] = []
            for BV in BVtry_lst:
                try:
                    teff = color_to_Teff('B-V',BV,FeH=feh,ref=ref,startype='dwarf')
                    BV_lst[ref].append(BV)
                    teff_lst[ref].append(teff)
                except:
                    continue

        ax = fig.add_axes([ifeh*0.24+0.05, 0.12, 0.218, 0.83])
        for ref in ref_lst:
            ax.plot(BV_lst[ref], teff_lst[ref], '-', label=ref)
        leg = ax.legend()
        for t in leg.get_texts():
            t.set_fontsize(8)
        leg.get_frame().set_alpha(0.1)
        x1, x2 = ax.get_xlim()
        y1, y2 = ax.get_ylim()
        ax.text(0.9*x1+0.1*x2, 0.9*y1+0.1*y2, '[Fe/H] = %3.1f'%feh)
        ax.set_xlabel('$B-V$')
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(10)
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(10)
        if ifeh == 0:
            ax.set_ylabel('$T_{\\rm eff}\,(\mathrm{K})$')
        else:
            ax.set_yticklabels([])
        ax.grid(True)

    fig.savefig('compare_BV_dwarf.png')

if __name__=='__main__':
    main()
