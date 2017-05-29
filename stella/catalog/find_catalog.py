from __future__ import print_function
import os
import struct
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

    nbyte, nrow, ncol, pos, dtype, fmtfunc = get_bintable_info(filename)

    if hip<=0 or hip > 118218 or hip > nrow:
        raise ItemNotFound(hip)

    infile = open(filename)
    infile.seek(pos+(hip-1)*nbyte, 0)
    item = fmtfunc(infile.read(nbyte))
    infile.close()

    if item['HIP'] <= 0:
        raise ItemNotFound(hip)

    return item

def find_TYC(starname):

    filename = os.path.join(os.getenv('STELLA_DATA'), 'catalog/TYC.fits')
    if not os.path.exists(filename):
        raise FileNotFond(filename)

    # get tyc1, tyc2, and tyc3
    if starname[0:3]=='TYC':
        g = starname[3:].split('-')
        tyc1 = int(g[0])
        tyc2 = int(g[1])
        tyc3 = int(g[2])

    # get information of FITS table
    nbyte, nrow, ncol, pos, dtype, fmtfunc = get_bintable_info(filename)

    target = np.int32((tyc1<<18) + (tyc2<<4) + (tyc3<<1))

    infile = open(filename)
    i1, i2 = 0, nrow-1
    while(i2-i1 > 1):
        i3 = int((i1+i2)/2)
        infile.seek(pos + i3*nbyte, 0)
        key = struct.unpack('>i',infile.read(4))[0]
        if target < key:
            i2 = i3
        elif target > key:
            i1 = i3
        else:
            break
    infile.seek(-4,1)
    row = fmtfunc(infile.read(nbyte))
    
    # looking for possible companion
    if i3 != nrow - 1:
        key2 = struct.unpack('>i',infile.read(4))[0]
        if key2 == key + 1:
            infile.seek(-4,1)
            row2 = fmtfunc(infile.read(nbyte))
            print('Warning: There are more than 1 star matched')
            t1 = (key2 & 0b11111111111111000000000000000000)/2**18
            t2 = (key2 & 0b00000000000000111111111111110000)/2**4
            t3 = (key2 & 0b00000000000000000000000000001110)/2
            print(t1,t2,t3,row2)
    infile.close()
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
