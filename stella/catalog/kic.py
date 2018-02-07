import os
import numpy as np
import astropy.io.fits as fits
from ..utils.fitsio import get_bintable_info
from ..utils.asciitable import structitem_to_dict
from .base import _get_KIC_number

class _KIC(object):
    '''Class for Kepler Input Catalog (`V/133
    <http://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=V/133>`_, Kepler
    Mission Team, 2009).

    The data file used in this function is complied from the 10th version of
    KIC. It contains 13,161,029 records with consecutive KIC numbers. Proper
    motions are available for 12,944,973 objects, or 98% of the entire sample.
    Parallaxes are provided for 958 objects, and physical parameters
    (*T*:sub:`eff`, log\ *g*, log\ *Z* and *R*) are available for 2,106,821
    objects, or 16% of the entire sample.
    For more details, see :ref:`Kepler Input Catalog<catalog_kic>`.

    '''
    def __init__(self):
        self.catfile = os.path.join(os.getenv('STELLA_DATA'), 'catalog/KIC.fits')
        self._data_info = None

    def _get_data_info(self):
        '''Get information of FITS table.'''
        nbyte, nrow, ncol, pos, dtype, fmtfunc = get_bintable_info(self.catfile)
        self._data_info = {
                'nbyte'  : nbyte,
                'nrow'   : nrow,
                'ncol'   : ncol,
                'pos'    : pos,
                'dtype'  : dtype,
                'fmtfunc': fmtfunc
                }
        
    def find_object(self, name, output='dict'):
        '''
        Find records in Kepler Input Catalog.

        .. csv-table:: Descriptions of returned parameters
            :header: Key, Type, Unit, Description
            :widths: 30, 30, 30, 120

            KIC,    integer32, ,       KIC number
            RAdeg,  float64,   deg,    Right ascension (*α*) at J2000
            DEdeg,  float64,   deg,    Declination (*δ*) at J2000
            pmRA,   float32,   mas/yr, Proper motion in Right ascension with cos(*δ*) factor
            pmDE,   float32,   mas/yr, Proper motion in Declination
            Plx,    float32,   mas,    Parallax
            umag,   float32,   mag,    *u* magnitude in SDSS system
            gmag,   float32,   mag,    *g* magnitude in SDSS system
            rmag,   float32,   mag,    *r* magnitude in SDSS system
            imag,   float32,   mag,    *i* magnitude in SDSS system
            zmag,   float32,   mag,    *z* magnitude in SDSS system
            grmag,  float32,   mag,    Magnitude in GRed band
            d51mag, float32,   mag,    Magnitude in DDO-51 filter
            kepmag, float32,   mag,    Magnitude in Kepler band
            flag_g, integer16, ,       "Galaxy flag (0 for star, 1 for galaxy)"
            flag_v, integer16, ,       "Variable flag (0 for normal, 1 for variable)"
            cq,     string5,   ,       Origin of Kepelr magnitude
            fv,     integer16, ,       "0 for outside Kepler FOV, 1/2 for inside, 2 for Kepler target"
            Teff,   integer16, K,      Effective temperature
            logg,   float32,   dex,    Surface gravity
            FeH,    float32,   dex,    Metallicity
            EBV,    float32,   mag,    Color excess in *B* − *V*
            Av,     float32,   mag,    Extinction in *V* magnitude
            R,      float32,   Rsun,   Stellar radius

        Args:
            name (string or integer): Name or number of star.
            output (string): Type of output results. Either *"dict"* or
                *"dtype"* (:class:`numpy.dtype`).
        Returns:
            dict or :class:`numpy.dtype`: Record in catalogue.
        Raises:
            UnrecognizedName: Input name can not be recognized
        Examples:
            Find *K*:sub:`p` magnitude of Kepler-13 (KOI-13, KIC 9941662)

            .. code-block:: python
        
                >>> from stella.catalog import KIC
                >>> res = KIC.find_object('KIC 9941662')
                >>> rec['kepmag']
                9.958000183105469

        '''

        kic = _get_KIC_number(name)

        if self._data_info is None:
            self._get_data_info()

        nrow    = self._data_info['nrow']
        nbyte   = self._data_info['nbyte']
        pos     = self._data_info['pos']
        fmtfunc = self._data_info['fmtfunc']

        infile = open(self.catfile, 'rb')

        if kic > 0 and kic <= nrow:
            infile.seek(pos+(kic-1)*nbyte,0)
            item = fmtfunc(infile.read(nbyte))
        else:
            pass

        infile.close()

        if output == 'ndarray':
            return item
        elif output == 'dict':
            return structitem_to_dict(item)
        else:
            return None

KIC = _KIC()

def plot_histogram():
    '''Plot magnitude and *T*:sub:`eff` - log\ *g* histograms of Kepler Input
    Catalogue.

    Examples:

        .. code-block:: python

            >>> from stella.catalog.kic import plot_histogram
            >>> plot_histogram()

    '''
    import matplotlib.pyplot as plt
    import matplotlib.ticker as tck

    filename = os.path.join(os.getenv('STELLA_DATA'), 'catalog/KIC.fits')
    data = fits.getdata(filename)
    mask = np.isnan(data['kepmag'])
    kps = data['kepmag'][~mask]
    
    mask1 = np.isnan(data['Teff'])
    mask2 = np.isnan(data['logg'])
    mask3 = (~mask1)*(~mask2)
    data3 = data[mask3]
    
    family = 'Sans Serif'
    label_fontsize=15
    tick_fontsize=12
    bins = np.arange(0,26)
    fig1 = plt.figure(figsize=(15,6), dpi=150)
    ax1 = fig1.add_axes([0.07,0.12,0.39,0.83])
    ax2 = fig1.add_axes([0.52,0.12,0.39,0.83])
    ax3 = fig1.add_axes([0.92,0.12,0.015,0.83])
    ax1.hist(kps, bins=bins, color='#1166aa', lw=0.5)
    ax1.set_yscale('log')
    ax1.set_xlabel('$K_\mathrm{p}$ Magnitude', fontsize=label_fontsize,family=family)
    ax1.set_ylabel('$N$', fontsize=label_fontsize,family=family)
    
    h, yedge, xedge = np.histogram2d(data3['logg'], data3['Teff'],
                        bins=[np.arange(0,6.01,0.1), np.arange(3000,12001,100)])
    cax = ax2.imshow(np.log10(h), cmap='Blues', interpolation='none', aspect='auto')
    cbar = fig1.colorbar(cax,cax=ax3)
    xticks = np.arange(0, 90+1, 20)
    ax2.set_xticks(xticks)
    ax2.set_xticklabels(np.int32(xedge[xticks]))
    yticks = np.arange(0, 60+1, 10)
    ax2.set_yticks(yticks)
    #for t in yticks:
    #    print(t, yedge[t])
    ax2.set_yticklabels(['%3.1f'%yedge[t] for t in yticks])
    ax2.set_xlim(90,0)
    ax2.set_ylim(60,0)
    ax2.set_xlabel('$T_\mathrm{eff}$ (K)',fontsize=label_fontsize,family=family)
    ax2.set_ylabel('$\log{g}$',fontsize=label_fontsize,family=family)
    cbar.set_label('$\log{N}$',fontsize=label_fontsize,family=family)
    #ax2.xaxis.set_major_locator(tck.MultipleLocator(1000))
    ax2.xaxis.set_minor_locator(tck.MultipleLocator(5))
    #ax2.yaxis.set_major_locator(tck.MultipleLocator(1.0))
    ax2.yaxis.set_minor_locator(tck.MultipleLocator(5))
    
    for ax in [ax1, ax2]:
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(tick_fontsize)
            tick.label1.set_family(family)
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(tick_fontsize)
            tick.label1.set_family(family)
    for tick in cbar.ax.get_yaxis().get_major_ticks():
        tick.label2.set_fontsize(tick_fontsize)
        tick.label2.set_family(family)
    fig1.savefig('catalogue_stat_KIC.png')

def convert_to_fits(path):
    '''Convert ASCII catalog files to FITS table.

    Args:
        path (string): Path to the KIC ASCII files.
    Returns:
        No returns.

    Examples:

        .. code-block:: python

            >>> import os
            >>> from stella.catalog.kic import convert_to_fits
            >>> astrodata = os.getenv('ASTRO_DATA')
            >>> path = os.path.join(astrodata, '/catalog/V/133/kic/')

    '''
    import time
    from .base import _str_to_float, _str_to_int

    filename_lst = []
    for direct in 'sn':
        for i in range(90):
            fn = os.path.join(path, '%s%02d.dat'%(direct,i))
            if os.path.exists(fn):
                filename_lst.append(fn)
    if len(filename_lst)==0:
        print('Error: Cannot find catalog file in %s'%path)

    #initialize types
    types = [
             ('KIC',    np.int32),
             ('RAdeg',  np.float64),
             ('DEdeg',  np.float64),
             ('pmRA',   np.float32),
             ('pmDE',   np.float32),
             ('Plx',    np.float32),
             ('umag',   np.float32),
             ('gmag',   np.float32),
             ('rmag',   np.float32),
             ('imag',   np.float32),
             ('zmag',   np.float32),
             ('grmag',  np.float32),
             ('d51mag', np.float32),
             #('Jmag',  np.float32),
             #('Hmag',  np.float32),
             #('Kmag',  np.float32),
             ('kepmag', np.float32),
             ('flag_g', np.int16),
             ('flag_v', np.int16),
             ('cq',     'S5'),
             ('fv',     np.int16),
             ('Teff',   np.int16),
             ('logg',   np.float32),
             ('FeH',    np.float32),
             ('EBV',    np.float32),
             ('Av',     np.float32),
             ('R',      np.float32),
            ]
    tmp = list(zip(*types))
    names = tmp[0]
    formats = tmp[1]
    record = np.dtype({'names':names,'formats':formats})

    data = []

    for ifile, fname in enumerate(filename_lst):
        t1 = time.time()
        data1 = []
        infile = open(fname)
        for row in infile:
            if row[0] == '#':
                continue
            kic    = int(row[0:8])
            ra     = np.float64(row[9:19])
            dec    = np.float64(row[19:29])
            pmra   = _str_to_float(row[30:38], np.NaN)
            pmde   = _str_to_float(row[39:47], np.NaN)
            plx    = _str_to_float(row[48:56], np.NaN)
            umag   = _str_to_float(row[57:63], np.NaN)
            gmag   = _str_to_float(row[63:70], np.NaN)
            rmag   = _str_to_float(row[70:77], np.NaN)
            imag   = _str_to_float(row[77:84], np.NaN)
            zmag   = _str_to_float(row[84:91], np.NaN)
            grmag  = _str_to_float(row[92:98], np.NaN)
            d51mag = _str_to_float(row[99:105], np.NaN)
            #magJ   = _str_to_float(row[106:112], np.NaN)
            #magH   = _str_to_float(row[112:119], np.NaN)
            #magK   = _str_to_float(row[119:126], np.NaN)
            kepmag = _str_to_float(row[127:133], np.NaN)
            flag_g = int(row[209])
            flag_v = int(row[211])
            cq     = row[213:218]
            fv     = int(row[223])
            teff   = _str_to_int(row[225:231], -1)
            logg   = _str_to_float(row[231:238], np.NaN)
            feh    = _str_to_float(row[238:245], np.NaN)
            EBV    = _str_to_float(row[246:252], np.NaN)
            Av     = _str_to_float(row[252:259], np.NaN)
            R      = _str_to_float(row[260:266], np.NaN)

            item = np.array((kic, ra, dec, pmra, pmde, plx,
                   umag, gmag, rmag, imag, zmag, grmag, d51mag, kepmag,
                   flag_g, flag_v, cq, fv,
                   teff, logg, feh, EBV, Av, R,
                   ),dtype=record)
            data1.append(item)
        infile.close()

        data1 = np.array(data1,dtype=record)

        data1 = np.sort(data1, order='KIC')

        t2 = time.time()
        dt = t2 - t1
        kic1 = data1[0]['KIC']
        kic2 = data1[-1]['KIC']
        print('%s %10d %10d %10d %10d %10.3f'%(os.path.basename(fname),
                kic1, kic2, kic2-kic1+1, data1.size, dt))

        for row in data1:
            data.append(row)

    data = np.array(data, dtype=record)

    pri_hdu = fits.PrimaryHDU()
    tbl_hdu = fits.BinTableHDU(data)
    hdu_lst = fits.HDUList([pri_hdu,tbl_hdu])
    outname = 'KIC.fits'
    if os.path.exists(outname):
        os.remove(outname)
    hdu_lst.writeto(outname)
