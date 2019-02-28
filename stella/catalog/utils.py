import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from astropy.coordinates import SkyCoord

def plot_skymap(ra, dec, figfile, projection='hammer', figsize=(8,4.5), dpi=150,
    size=1, alpha=0.1):
    """Plot a skymap of sample stars.

    Args:
        ra (list or :class:`numpy.ndarray`): List of RA.
        dec (list or :class:`numpy.ndarray`): List of Dec.
        figfile (str): Name of output figure.
        projection (str): Projection of map. Avialable options include
            'aitoff', 'mollweide', and 'hammer'.
        figsize (tuple): Size of figure in tuple (width, height).
        dpi (int): DPI of figure.
        size (float): Size of data points.
        alpha (float): A float between (0, 1] representing the transparency of
            data points.
    """

    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax  = fig.add_axes([0.05, 0.05, 0.9, 0.9], projection=projection)
    ra  = np.deg2rad(180-ra)
    dec = np.deg2rad(dec)
    ax.scatter(ra, dec, s=size, c='#1166aa', alpha=alpha, lw=0)
    ax.grid(True, color='k', linestyle=':', linewidth=0.5)
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

def plot_histogram(x, xlabel, figfile, bins, figsize=(8,6), dpi=150,
    color='#1166aa', alpha=1, rwidth=0.9, yscale='log', xlim=None, ylim=None,
    xticks=None, ticksize=13, labelsize=15):
    """Plot a histogram.

    Args:
        x (list or :class:`numpy.ndarray`): List of data.
        xlabel (str): Label in X-axis.
        figfile (str): Name of output figure.
        bins (list or :class:`numpy.ndarray`): Bins of data.
        figsize (tuple): Size of figure in tuple (width, height).
        dpi (int): DPI of figure.
        color (str): Color of histogram bars.
        alpha (float): A float between (0, 1] representing the transparency of
            histogram bars.
        rwidth (float): Relative width of histogram bars.
        yscale (str): Scale of Y axis. Either 'linear' or 'log'.
        xlim (tuple): Limits of X axis.
        ylim (tuple): Limits of Y axis.
        xticks (list): List of ticks in X axis.
        ticksize (int): Size of tick labels.
        labelsize (int): Size of X and Y axis labels.
    """

    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax = fig.add_axes([0.1,0.1,0.88,0.85])
    ax.hist(x, bins=bins, color=color, alpha=alpha, rwidth=rwidth)
    ax.set_axisbelow(True)
    ax.set_facecolor('#dddddd')
    ax.yaxis.grid(True, color='w', linestyle='-', linewidth=1)
    ax.set_yscale(yscale)

    # change tick size
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(ticksize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(ticksize)

    if xlim is None:
        xlim = bins[0], bins[-1]
    ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    if xticks is not None:
        ax.set_xticks(xticks)
    ax.set_xlabel(xlabel, fontsize=labelsize)
    ax.set_ylabel('$N$', fontsize=labelsize)

    # save the figure
    fig.savefig(figfile)
    plt.close(fig)

def plot_histogram2d(x, y, xbins, ybins, xlabel, ylabel, figfile,
    figsize=(8,6), dpi=150, reverse_x=False, reverse_y=False, scale='log',
    ticksize=13, labelsize=15):
    """Plot a 2-D histogram of H-R diagram.

    Args:
        x (list or :class:`numpy.ndarray`): List of x data.
        y (list or :class:`numpy.ndarray`): List of y data.
        xbins (list or :class:`numpy.ndarray`): Bins of data along x axis.
        ybins (list or :class:`numpy.ndarray`): Bins of data along y axis.
        xlabel (str): Label in x-axis.
        ylabel (str): Label in y-axis.
        figfile (str): Name of output figure.
        figsize (tuple): Size of figure in tuple (width, height).
        dpi (int): DPI of figure.
        reverse_x (bool): Reverse x axis if *True*.
        reverse_y (bool): Reverse y axis if *True*.
        scale (str): Scale of y axis.
        ticksize (int): Size of tick labels.
        labelsize (int): Size of X and Y axis labels.
    """

    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax1 = fig.add_axes([0.1,0.1,0.75,0.85])
    ax2 = fig.add_axes([0.88,0.1,0.03,0.85])

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

    ax1.set_xlabel(xlabel, fontsize=labelsize)
    ax1.set_ylabel(ylabel, fontsize=labelsize)
    cbar.set_label('$N$', fontsize=labelsize)

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

    # change tick size
    for tick in ax1.xaxis.get_major_ticks():
        tick.label1.set_fontsize(ticksize)
    for tick in ax1.yaxis.get_major_ticks():
        tick.label1.set_fontsize(ticksize)
    cbar.ax.tick_params(labelsize=ticksize)

    fig.savefig(figfile)
    plt.close(fig)
