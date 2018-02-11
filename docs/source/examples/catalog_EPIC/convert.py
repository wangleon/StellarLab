#!/usr/bin/env python3
import os
import astropy.io.fits as fits
import numpy as np

inputfile = '/opt/ref/data/ApJS/224/2/table5.dat'

types = [
        ('EPIC',     np.int32),
        ('Teff',     np.int16),
        ('E_Teff',   np.int16),
        ('e_Teff',   np.int16),
        ('logg',     np.float32),
        ('E_logg',   np.float32),
        ('e_logg',   np.float32),
        ('FeH',      np.float32),
        ('E_FeH',    np.float32),
        ('e_FeH',    np.float32),
        ('Rad',      np.float32),
        ('E_Rad',    np.float32),
        ('e_Rad',    np.float32),
        ('Mass',     np.float32),
        ('E_Mass',   np.float32),
        ('e_Mass',   np.float32),
        ('rho',      np.float32),
        ('E_rho',    np.float32),
        ('e_rho',    np.float32),
        ('Dist',     np.float32),
        ('E_Dist',   np.float32),
        ('e_Dist',   np.float32),
        ('E(B-V)',   np.float32),
        ('E_E(B-V)', np.float32),
        ('e_E(B-V)', np.float32),
        ('Flag',     'S3'),
        ('RAdeg',    np.float64),
        ('DEdeg',    np.float64),
        ]

tmp = list(zip(*types))
record = np.dtype({'names':tmp[0],'formats':tmp[1]})

data = []
prev_epic = 0
infile = open(inputfile)
for row in infile:
    epic = int(row[0:9])
    teff = int(row[10:15])
    E_teff = int(row[16:20])
    e_teff = int(row[21:25])
    logg   = float(row[26:32])
    E_logg = float(row[33:38])
    e_logg = float(row[39:44])
    feh    = float(row[45:51])
    E_feh  = float(row[52:57])
    e_feh  = float(row[58:63])
    rad    = float(row[64:71])
    E_rad  = float(row[72:79])
    e_rad  = float(row[80:87])
    mass   = float(row[88:93])
    E_mass = float(row[94:99])
    e_mass = float(row[100:105])
    rho    = float(row[106:115])
    E_rho  = float(row[116:125])
    e_rho  = float(row[126:135])
    dist   = float(row[136:145])
    E_dist = float(row[146:155])
    e_dist = float(row[156:165])
    ebv    = float(row[166:173])
    E_ebv  = float(row[174:180])
    e_ebv  = float(row[181:187])
    flag   = row[188:191]
    radeg  = float(row[192:203])
    dedeg  = float(row[204:215])

    item = np.array((epic, teff, E_teff, e_teff, logg, E_logg, e_logg,
                           feh,  E_feh,  e_feh,  rad,  E_rad,  e_rad,
                           mass, E_mass, e_mass, rho,  E_rho,  e_rho,
                           dist, E_dist, e_dist, ebv,  E_ebv,  e_ebv,
                           flag, radeg, dedeg),dtype=record)

    data.append(item)

    if epic <= prev_epic:
        print epic
    prev_epic = epic

infile.close()

data = np.array(data, dtype=record)

pri_hdu = fits.PrimaryHDU()
tbl_hdu = fits.BinTableHDU(data)
hdu_lst = fits.HDUList([pri_hdu,tbl_hdu])
outputfile = 'EPIC.fits'
if os.path.exists(outputfile):
    os.remove(outputfile)
hdu_lst.writeto(outputfile)
