#!/usr/bin/env python3
import os
import numpy as np
import astropy.io.fits as fits
from stella.catalog.base import _str_to_float

inputfile = os.path.join(os.getenv('ASTRO_DATA'), 'catalog/I/239/hip_main.dat')

types = [
        ('HIP',    np.int32),
        ('RAdeg',  np.float64),
        ('DEdeg',  np.float64),
        ('Vmag',   np.float32),
        ('Plx',    np.float32),
        ('e_Plx',  np.float32),
        ('pmRA',   np.float32),
        ('pmDE',   np.float32),
        ('e_pmRA', np.float32),
        ('e_pmDE', np.float32),
        ('BTmag',  np.float32),
        ('e_BTmag',np.float32),
        ('VTmag',  np.float32),
        ('e_VTmag',np.float32),
        ('B-V',    np.float32),
        ('e_B-V',  np.float32),
        ('r_B-V',  'S1'),
        ('V-I',    np.float32),
        ('e_V-I',  np.float32),
        ('r_V-I',  'S1'),
        ('Hpmag',  np.float32),
        ('e_Hpmag',np.float32),
        ('Hpscat', np.float32),
        ('o_Hpmag',np.int16),
        #('CCDM',   'S10'),
        #('HD',     np.int32),
        #('BD',     'S10'),
        #('CoD',    'S10'),
        ('SpType', 'S12'),
        ('r_SpType','S1'),
        ]
tmp = list(zip(*types))
record = np.dtype({'names':tmp[0],'formats':tmp[1]})

fill_item = np.array((0, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN,
                         np.NaN, np.NaN, np.NaN, np.NaN, np.NaN,
                         np.NaN, np.NaN, np.NaN, np.NaN, np.NaN,
                         '',     np.NaN, np.NaN, '', np.NaN, np.NaN, np.NaN,
                         -32768, '',''),dtype=record)

data = {}

infile = open(inputfile)
for row in infile:
    hip = int(row[8:14])
    rah = int(row[17:19])
    ram = int(row[20:22])
    ras = float(row[23:28])
    radeg = (rah + ram/60. + ras/3600.)*15.
    ded = abs(int(row[30:32]))
    dem = int(row[33:35])
    des = float(row[36:40])
    dedeg = ded + dem/60. + des/3600.
    if row[29]=='-':
        dedeg = -dedeg
    vmag     = _str_to_float(row[41:46], np.NaN)
    plx      = _str_to_float(row[79:86], np.NaN)
    e_plx    = _str_to_float(row[119:125], np.NaN)
    pmRA     = _str_to_float(row[87:95], np.NaN)
    pmDE     = _str_to_float(row[96:104], np.NaN)
    e_pmRA   = _str_to_float(row[126:132], np.NaN)
    e_pmDE   = _str_to_float(row[133:139], np.NaN)
    BTmag    = _str_to_float(row[217:223], np.NaN)
    VTmag    = _str_to_float(row[230:236], np.NaN)
    e_BTmag  = _str_to_float(row[224:229], np.NaN)
    e_VTmag  = _str_to_float(row[237:242], np.NaN)
    BV       = _str_to_float(row[245:251], np.NaN)
    e_BV     = _str_to_float(row[252:257], np.NaN)
    r_BV     = row[258].strip()
    VI       = _str_to_float(row[260:264], np.NaN)
    e_VI     = _str_to_float(row[265:269], np.NaN)
    r_VI     = row[270].strip()
    Hpmag    = _str_to_float(row[274:281], np.NaN)
    e_Hpmag  = _str_to_float(row[282:288], np.NaN)
    Hpscat   = _str_to_float(row[289:294], np.NaN)
    if row[295:298].strip()=='':
        o_Hpmag = 0
    else:
        o_Hpmag = int(row[295:298])

    if not np.isnan(pmRA):
        pm_ra = pmRA*1e-3/3600. # convert pm_RA  from mas/yr to degree/yr
        #radeg += (2000.0-1991.25)*pm_ra/math.cos(dedeg/180.*math.pi)

    if not np.isnan(pmDE):
        pm_de = pmDE*1e-3/3600. # convert pm_Dec from mas/yr to degree/yr
        #dedeg += (2000.0-1991.25)*pm_de

    SpType   = row[435:447].strip()
    r_SpType = row[448].strip()

    item = np.array((hip, radeg, dedeg, vmag, plx, e_plx,
                     pmRA, pmDE, e_pmRA, e_pmDE,
                     BTmag, e_BTmag, VTmag, e_VTmag,
                     BV, e_BV, r_BV, VI, e_VI, r_VI,
                     Hpmag, e_Hpmag, Hpscat, o_Hpmag,
                     #CCDM, HD, BD, CoD,
                     SpType, r_SpType,
                    ),dtype=record)
    if hip in data:
        print('Error: Duplicate Records for HIP', hip)
    data[hip] = item
infile.close()

newdata = []
for hip in range(1, max(data.keys())+1):
    if hip in data:
        newdata.append(data[hip])
    else:
        newdata.append(fill_item)
newdata = np.array(newdata, dtype=record)

pri_hdu = fits.PrimaryHDU()
tbl_hdu = fits.BinTableHDU(newdata)
hdu_lst = fits.HDUList([pri_hdu,tbl_hdu])

outputfile='HIP.fits'
if os.path.exists(outputfile):
    os.remove(outputfile)
hdu_lst.writeto(outputfile)
