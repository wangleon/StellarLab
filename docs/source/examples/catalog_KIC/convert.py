#!/usr/bin/env python3
import os
import time
import numpy as np
import astropy.io.fits as fits
from stella.catalog.base import _str_to_float, _str_to_int

path = os.path.join(os.getenv('ASTRO_DATA'), 'catalog/V/133/kic/')

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

outputfile = 'KIC.fits'
if os.path.exists(outputfile):
    os.remove(outputfile)
hdu_lst.writeto(outputfile)
