#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from stella.parameter.teff import color_to_Teff

def main():

    BVtry_lst = np.arange(0.0,1.8+1e-4,0.01)
    ref_lst = ['Alonso1996','Ramirez2005','GB2009','Casagrande2010']

    fig = plt.figure(figsize=(8,8),dpi=150,tight_layout=True)
    for ifeh, feh in enumerate([-2,-1,0,0.2]):
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

        ax = fig.add_subplot(221+ifeh)
        for ref in ref_lst:
            ax.plot(BV_lst[ref], teff_lst[ref], '-', label=ref)
        leg = ax.legend()
        for t in leg.get_texts():
            t.set_fontsize(8)
        leg.get_frame().set_alpha(0.1)
        ax.set_xlabel('$B-V$')
        ax.set_ylabel('$T_{\\rm eff}$')
        x1, x2 = ax.get_xlim()
        y1, y2 = ax.get_ylim()
        ax.text(0.9*x1+0.1*x2, 0.9*y1+0.1*y2, '[Fe/H] = %3.1f'%feh)
        ax.grid(True)

    fig.savefig('compare_BV.png')
    plt.show()

if __name__=='__main__':
    main()
