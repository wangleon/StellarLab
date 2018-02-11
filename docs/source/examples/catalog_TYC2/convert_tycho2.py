#!/usr/bin/env python3
import os
import numpy as np
import astropy.io.fits as fits
from stella.catalog.base import _str_to_float, _str_to_int

path = os.path.join(os.getenv('ASTRO_DATA'), 'catalog/I/259/')

# load duplicate list
dup_tyc_lst = []
infile1 = open('duplicate.txt')
for row in infile1:
    g = row.split('|')[0].split()
    tyc1 = int(g[0])
    tyc2 = int(g[1])
    tyc3 = int(g[2])
    dup_tyc_lst.append((tyc1,tyc2,tyc3))
infile1.close()

types = [
        ('TYC',      np.int32),
        ('RAdeg',    np.float64),
        ('DEdeg',    np.float64),
        ('pmRA',     np.float32),
        ('pmDE',     np.float32),
        ('e_RA',     np.int16),
        ('e_DE',     np.int16),
        ('e_pmRA',   np.float32),
        ('e_pmDE',   np.float32),
        ('epoch_RA', np.float32),
        ('epoch_DE', np.float32),
        ('num',      np.int16),
        ('q_RA',     np.float32),
        ('q_DE',     np.float32),
        ('q_pmRA',   np.float32),
        ('q_pmDE',   np.float32),
        ('BTmag',    np.float32),
        ('e_BTmag',  np.float32),
        ('VTmag',    np.float32),
        ('e_VTmag',  np.float32),
        ('prox',     np.float32),
        ]
tmp = zip(*types)
record = np.dtype({'names':tmp[0],'formats':tmp[1]})

fn_lst = ['tyc2.dat.%02d'%i for i in range(20)]
data = []
for fn in fn_lst:
    print fn
    infile = open(os.path.join(path, fn))
    for row in infile:
        tyc1   = int(row[0:4])
        tyc2   = int(row[5:10])
        tyc3   = int(row[11])
        tyc = np.int32((tyc1<<18) + (tyc2<<4) + (tyc3<<1))

        radeg    = _str_to_float(row[15:27], np.Nan)
        dedeg    = _str_to_float(row[28:40], np.Nan)
        pmRA     = _str_to_float(row[41:48], np.Nan)
        pmDE     = _str_to_float(row[49:56], np.Nan)
        e_ra     = _str_to_int(row[57:60], -1)
        e_de     = _str_to_int(row[61:64], -1)
        e_pmRA   = _str_to_float(row[65:69], np.Nan)
        e_pmDE   = _str_to_float(row[70:74], np.Nan)
        epoch_ra = _str_to_float(row[75:82], np.Nan)
        epoch_de = _str_to_float(row[83:90], np.Nan)
        num      = _str_to_int(row[91:93], -1)
        q_RA     = _str_to_float(row[94:97], np.Nan)
        q_DE     = _str_to_float(row[98:101], np.Nan)
        q_pmRA   = _str_to_float(row[102:105], np.Nan)
        q_pmDE   = _str_to_float(row[106:109], np.Nan)
        BTmag    = _str_to_float(row[110:116], np.Nan)
        e_BTmag  = _str_to_float(row[117:122], np.Nan)
        VTmag    = _str_to_float(row[123:129], np.Nan)
        e_VTmag  = _str_to_float(row[130:135], np.Nan)
        prox     = _str_to_float(row[136:139], np.Nan)*0.1

        item = np.array((tyc, radeg, dedeg, pmRA, pmDE, 
                         e_ra, e_de, e_pmRA, e_pmDE,
                         epoch_ra, epoch_de, num,
                         q_RA, q_DE, q_pmRA, q_pmDE,
                         BTmag, e_BTmag, VTmag, e_VTmag, prox), dtype=record)
        data.append(item)
    infile.close()


for fn in ['suppl_1.dat', 'suppl_2.dat']:
    infile = open(os.path.join(path, fn))
    for row in infile:
        tyc1   = int(row[0:4])
        tyc2   = int(row[5:10])
        tyc3   = int(row[11])

        if (tyc1,tyc2,tyc3) in dup_tyc_lst:
            m = 1
        else:
            m = 0

        tyc = np.int32((tyc1<<18) + (tyc2<<4) + (tyc3<<1) + m)

        radeg    = _str_to_float(row[15:27], np.Nan)
        dedeg    = _str_to_float(row[28:40], np.Nan)
        pmRA     = _str_to_float(row[41:48], np.Nan)
        pmDE     = _str_to_float(row[49:56], np.Nan)
        e_ra     = _str_to_int(row[57:60], -1)
        e_de     = _str_to_int(row[61:64], -1)
        e_pmRA   = _str_to_float(row[65:69], np.Nan)
        e_pmDE   = _str_to_float(row[70:74], np.Nan)
        epoch_ra = _str_to_float(row[75:82], np.Nan)
        epoch_de = _str_to_float(row[83:90], np.Nan)
        num      = _str_to_int(row[91:93], -1)
        q_RA     = _str_to_float(row[94:97], np.Nan)
        q_DE     = _str_to_float(row[98:101], np.Nan)
        q_pmRA   = _str_to_float(row[102:105], np.Nan)
        q_pmDE   = _str_to_float(row[106:109], np.Nan)
        BTmag    = _str_to_float(row[110:116], np.Nan)
        e_BTmag  = _str_to_float(row[117:122], np.Nan)
        VTmag    = _str_to_float(row[123:129], np.Nan)
        e_VTmag  = _str_to_float(row[130:135], np.Nan)
        prox     = _str_to_float(row[136:139], np.Nan)*0.1

        item = np.array((tyc, radeg, dedeg, pmRA, pmDE,
                         e_ra, e_de, e_pmRA, e_pmDE,
                         epoch_ra, epoch_de, num,
                         q_RA, q_DE, q_pmRA, q_pmDE,
                         BTmag, e_BTmag, VTmag, e_VTmag, prox), dtype=record)
        data.append(item)
    infile.close()

data = np.array(data,dtype=record)
data = np.sort(data, order='TYC')

pri_hdu = fits.PrimaryHDU()
tbl_hdu = fits.BinTableHDU(data)
hdu_lst = fits.HDUList([pri_hdu,tbl_hdu])
outputfile = 'TYC2.fits'
if os.path.exists(outputfile):
    os.remove(outputfile)
hdu_lst.writeto(outputfile)

data = fits.getdata(outputfile)
tyc_lst = data['TYC']
diff = np.diff(tyc_lst)
mask = diff==0
for i, m in enumerate(mask):
    if m:
        print(data[i]['TYC'], data[i]['Note'], data[i+1]['TYC'],data[i+1]['Note'])

