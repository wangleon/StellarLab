import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from astropy.coordinates import SkyCoord

def plot_skymap(ra, dec, figfile, projection='hammer', figsize=(8,4.5), dpi=150,
    size=1, alpha=0.1):
    '''Plot a skymap of sample stars.

    Args:
        ra (list or :class:`numpy.array`): List of RA.
        dec (list or :class:`numpy.array`): List of Dec.
        figfile (string): Name of output figure.
        projection (string): Projection of map. Avialable options include
            'aitoff', 'mollweide', and 'hammer'.
        figsize (tuple): Size of figure in tuple (width, height).
        dpi (integer): DPI of figure.
        size (float): Size of data points.
        alpha (float): A float between (0, 1] representing the transparency of
            data points.
    Returns:
        No returns.
    '''

    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax  = fig.add_axes([0.05, 0.05, 0.9, 0.9], projection=projection)
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
    c = SkyCoord(np.arange(360), 0, frame='galactic', unit='deg')
    eq = c.transform_to('icrs')
    ra_lst  = eq.ra.degree
    dec_lst = eq.dec.degree
    imaxdiff = np.abs(np.diff(ra_lst)).argmax()
    ra_lst  = np.roll(ra_lst, -imaxdiff-1)
    dec_lst = np.roll(dec_lst, -imaxdiff-1)
    ra_lst  = np.deg2rad(180-ra_lst)
    dec_lst = np.deg2rad(dec_lst)
    ax.plot(ra_lst, dec_lst, 'k--', lw=0.8, label='Galactic Plane')

    ax.legend(loc='upper right')

    # save the figure
    fig.savefig(figfile)
    plt.close(fig)

def plot_histogram(mags, xlabel, figfile, bins, figsize=(8,6), dpi=150,
    color='#1166aa', alpha=1, linewidth=0.5, yscale='log'):
    '''Plot a histogram.

    Args:
        mags (list or :class:`numpy.array`): List of magnitudes.
        xlabel (string): Label in X-axis.
        figfile (string): Name of output figure.
        bins (list or :class:`numpy.array`): Bins of magnitudes.
        figsize (tuple): Size of figure in tuple (width, height).
        dpi (integer): DPI of figure.
        color (string): Color of histogram bars.
        alpha (float): A float between (0, 1] representing the transparency of
            histogram bars.
        linewidth (float): Line width of histogram bars.
        yscale (string): Scale of Y axis. Either 'linear' or 'log'.
    Returns:
        No returns.
    '''

    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax = fig.add_axes([0.1,0.1,0.85,0.85])
    ax.hist(mags, bins=bins, color=color, alpha=alpha, lw=linewidth)
    ax.set_yscale(yscale)
    ax.set_xlabel(xlabel)
    ax.set_ylabel('$N$')

    # save the figure
    fig.savefig(figfile)
    plt.close(fig)

def plot_histogram2d(x, y, x_range, y_range, xlabel, ylabel, figfile,
    reverse_x = False, reverse_y = False, figsize=(8,6), dpi=150, scale='log'):
    '''Plot a 2-D histogram of H-R diagram.

    Args:
        x (list or :class:`numpy.array`): List of x data.
        y (list or :class:`numpy.array`): List of y data.
        x_range (tuple): Beginning, end, and step of x data: (x1, x2, dx).
        y_range (tuple): Beginning, end, and step of y data: (y1, y2, dy).
        xlabel (string): Label in x-axis.
        ylabel (string): Label in y-axis.
        figfile (string): Name of output figure.
        figsize (tuple): Size of figure in tuple (width, height).
        dpi (integer): DPI of figure.
    Returns:
        No returns.
    '''

    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax1 = fig.add_axes([0.1,0.1,0.75,0.85])
    ax2 = fig.add_axes([0.88,0.1,0.03,0.85])

    x1, x2, dx = x_range
    y1, y2, dy = y_range
    xbins = np.arange(x1, x2+1e-6, dx)
    ybins = np.arange(y1, y2+1e-6, dy)
    nx = xbins - 1
    ny = ybins - 1
    
    if scale=='log':
        norm = mcolors.LogNorm()
    else:
        norm = mcolors.Normalize()

    _,_,_,cax = ax1.hist2d(x, y, bins=(xbins, ybins), cmap='Blues', norm=norm)
    cbar = fig.colorbar(cax, cax=ax2)

    if reverse_x:
        _x1, _x2 = ax1.get_xlim()
        ax1.set_xlim(_x2, _x1)
    if reverse_y:
        _y1, _y2 = ax1.get_ylim()
        ax1.set_ylim(_y2, _y1)

    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    cbar.set_label('$N$')

    if scale=='log':
        c1, c2 = cbar.get_clim()
        mticks = []
        power1 = math.floor(math.log10(c1))
        power2 = math.ceil(math.log10(c2))
        for power in np.arange(power1, power2):
            for v in np.arange(1,10)*10**power:
                if c1 <= v <= c2:
                    mticks.append(v)
        mticks = cax.norm(mticks)
        cbar.ax.yaxis.set_ticks(mticks, minor=True)
    else:
        cbar.ax.minorticks_on()

    fig.savefig(figfile)
    plt.close(fig)
