import os
import numpy as np
import astropy.io.fits as fits

from ..utils.fitsio import get_bintable_info
from .errors import FileNotFound, ItemNotFound

def find_HIP(starname):

    filename = os.path.join(os.getenv('STELLA_DATA'), 'catalog/HIP.fits')
    if not os.path.exists(filename):
        raise FileNotFond(filename)

    try:
        hip = int(starname)
    except:
        if starname[0:3] == 'HIP':
            hip = int(starname[3:])

    if hip <= 0:
        raise ItemNotFound(hip)

    nbyte, nrow, ncol, pos, dtype, fmtfunc = get_bintable_info(filename)

    if hip > nrow:
        raise ItemNotFound(hip)

    infile = open(filename)
    infile.seek(pos+(hip-1)*nbyte,0)
    item = fmtfunc(infile.read(nbyte))
    infile.close()

    if item['HIP'] <=0:
        raise ItemNotFound(hip)

    return item

def find_HIP2(starname):


    filename = os.path.join(os.getenv('STELLA_DATA'), 'catalog/HIP2.fits')
    if not os.path.exists(filename):
        raise FileNotFond(filename)

    if type(starname)==type(1):
        hip = starname
    elif starname[0:3]=='HIP':
        hip = int(starname[3:])
    elif starname.isdigit():
        hip = int(starname)

    if hip<=0 or hip > 118218:
        return None

    f = fits.open(filename)
    data = f[1].data
    i = np.searchsorted(data['HIP'],hip)
    row = data[i]
    f.close()
    return row

def find_TYC(starname):

    filename = os.path.join(os.getenv('STELLA_DATA'), 'catalog/TYC.fits')
    if not os.path.exists(filename):
        raise FileNotFond(filename)

    if starname[0:3]=='TYC':
        g = starname[3:].split('-')
        tyc1 = int(g[0])
        tyc2 = int(g[1])
        tyc3 = int(g[2])


    f = fits.open(filename)
    data = f[1].data
    f.close()
    m1 = data['TYC1']==tyc1
    m2 = data[m1]['TYC2']==tyc2
    m3 = data[m1][m2]['TYC3']==tyc3
    row = data[m1][m2][m3][0]
    return row

def find_KIC(starname):

    filename = os.path.join(os.getenv('STELLA_DATA'), 'catalog/KIC.fits')
    if not os.path.exists(filename):
        raise FileNotFond(filename)

    try:
        kic = int(starname)
    except:
        if starname[0:3]=='KIC':
            kic = int(starname[3:])

    if kic <= 0:
        raise ItemNotFound

    nbyte, nrow, ncol, pos, dtype, fmtfunc = get_bintable_info(filename)

    if kic > nrow:
        raise ItemNotFound(kic)

    infile = open(filename)
    infile.seek(pos+(kic-1)*nbyte,0)
    item = fmtfunc(infile.read(nbyte))
    infile.close()

    return item
