#!/usr/bin/env python3
import os
import numpy as np
import astropy.io.fits as fits
from stella.catalog.base import _str_to_float

inputfile = os.path.join(os.getenv('ASTRO_DATA'), 'catalog/I/311/hip2.dat')

types = [
        ('HIP',    np.int32),
        ('RAdeg',  np.float64),
        ('DEdeg',  np.float64),
        ('Plx',    np.float32),
        ('e_Plx',  np.float32),
        ('pmRA',   np.float32),
        ('pmDE',   np.float32),
        ('e_pmRA', np.float32),
        ('e_pmDE', np.float32),
        ('B-V',    np.float32),
        ('e_B-V',  np.float32),
        ('V-I',    np.float32),
        ('Hpmag',  np.float32),
        ('e_Hpmag',np.float32),
        ('Hpscat', np.float32),
        ]
tmp = list(zip(*types))
record = np.dtype({'names':tmp[0],'formats':tmp[1]})

fill_item = np.array(
                (0, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN,
                    np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN),
                dtype=record)

data = {}

infile = open(inputfile)
for row in infile:
    hip     = int(row[0:6])
    radeg   = float(row[15:28])/math.pi*180.
    dedeg   = float(row[29:42])/math.pi*180.
    plx     = _str_to_float(row[43:50], np.NaN)
    e_plx   = _str_to_float(row[83:89], np.NaN)
    pmRA    = _str_to_float(row[51:59], np.NaN)
    pmDE    = _str_to_float(row[60:68], np.NaN)
    e_pmRA  = _str_to_float(row[90:96], np.NaN)
    e_pmDE  = _str_to_float(row[97:103], np.NaN)
    Hpmag   = _str_to_float(row[129:136], np.NaN)
    e_Hpmag = _str_to_float(row[137:143], np.NaN)
    Hpscat  = _str_to_float(row[144:149], np.NaN)
    BV      = _str_to_float(row[152:158], np.NaN)
    e_BV    = _str_to_float(row[159:164], np.NaN)
    VI      = _str_to_float(row[165:171], np.NaN)

    if not np.isnan(pmRA):
        pm_ra = pmRA*1e-3/3600. # convert pm_RA  from mas/yr to degree/yr
        #radeg += (2000.0-1991.25)*pm_ra/math.cos(dedeg/180.*math.pi)

    if not np.isnan(pmDE):
        pm_de = pmDE*1e-3/3600. # convert pm_Dec from mas/yr to degree/yr
        #dedeg += (2000.0-1991.25)*pm_de

    item = np.array((hip,radeg, dedeg, plx, e_plx,
                    pmRA, pmDE, e_pmRA, e_pmDE,
                    BV, e_BV, VI, Hpmag,e_Hpmag, Hpscat), dtype=record)
    if hip in data:
        print('Error: Duplicate record for HIP', hip)
    data[hip] = item

infile.close()

count_good, count_null = 0, 0
newdata = []
for hip in range(1, max(data.keys())+1):
    if hip in data:
        newdata.append(data[hip])
        count_good += 1
    else:
        print(hip)
        newdata.append(fill_item)
        count_null += 1
newdata = np.array(newdata, dtype=record)

pri_hdu = fits.PrimaryHDU()
tbl_hdu = fits.BinTableHDU(newdata)
hdu_lst = fits.HDUList([pri_hdu,tbl_hdu])

outputfile='HIP2.fits'
if os.path.exists(outputfile):
    os.remove(outputfile)
hdu_lst.writeto(outputfile)
