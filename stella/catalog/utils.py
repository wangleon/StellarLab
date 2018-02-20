import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord

def plot_skymap(ra, dec, figfile, figsize=(8,4.5), dpi=150, size=1, alpha=0.1):
    '''Plot a skymap of sample stars.

    Args:
        ra (list or :class`numpy.array)`: List of RA.
        dec (list or :class`numpy.array)`: List of Dec.
        figfile (string): Name of output figure.
        figsize (tuple): Size of figure in tuple (width, height).
        dpi (integer): DPI of figure.
        size (float): Size of data points.
        alpha (float): A float between (0, 1] representing the transparency of
            data points.
    Returns:
        No returns.
    '''

    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax  = fig.add_axes([0.05, 0.05, 0.9, 0.9], projection='hammer')
    ra  = np.deg2rad(180-ra)
    dec = np.deg2rad(dec)
    ax.scatter(ra, dec, s=size, c='#1166aa', alpha=alpha, lw=0)
    ax.grid(True)
    ax.set_xticklabels(['%dh'%v for v in np.arange(22,0,-2)])
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(10)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(10)

    # plot the galactic plane
    ra_lst = []
    dec_lst = []
    #c = SkyCoord(0,0, frame='icrs', unit='deg')
    for l in np.arange(0,360):
        c = SkyCoord(l, 0, frame='galactic', unit='deg')
        eq = c.transform_to('icrs')
        ra_lst.append(eq.ra.degree)
        dec_lst.append(eq.dec.degree)
    ra_lst = np.array(ra_lst)
    dec_lst = np.array(dec_lst)
    imaxdiff = np.abs(np.diff(ra_lst)).argmax()
    ra_lst = np.roll(ra_lst, -imaxdiff-1)
    dec_lst = np.roll(dec_lst, -imaxdiff-1)
    ra_lst = np.deg2rad(180-ra_lst)
    dec_lst = np.deg2rad(dec_lst)
    ax.plot(ra_lst, dec_lst, 'k--', lw=0.8)

    # save the figure
    fig.savefig(figfile)
    plt.close(fig)
