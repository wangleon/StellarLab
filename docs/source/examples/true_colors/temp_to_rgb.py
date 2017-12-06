#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
from stella.utils.vision import temp_to_rgb

def main():

    temp_lst = np.arange(0, 50000, 50)
    r_lst, g_lst, b_lst = [], [], []
    for temp in temp_lst:
        r,g,b = temp_to_rgb(temp)
        r_lst.append(r)
        g_lst.append(g)
        b_lst.append(b)
    
    fig = plt.figure(dpi=150, figsize=(8,4), tight_layout=True)
    ax = fig.gca()
    ax.plot(temp_lst, r_lst, 'r-')
    ax.plot(temp_lst, g_lst, 'g-')
    ax.plot(temp_lst, b_lst, 'b-')
    ax.set_ylim(-0.1, 1.1)
    ax.xaxis.set_major_locator(tck.MultipleLocator(5000))
    ax.xaxis.set_minor_locator(tck.MultipleLocator(1000))
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel('RGB Color')
    fig.savefig('temp_rgb.png')
    plt.show()

if __name__=='__main__':
    main()
