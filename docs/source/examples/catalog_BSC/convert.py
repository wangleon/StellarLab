#!/usr/bin/env python3
import os
import numpy as np
import astropy.io.fits as fits
from stella.catalog.base import _str_to_float, _str_to_int

def main():
    inputfile = os.path.join(os.getenv('ASTRO_DATA'), 'catalog/V/50/catalog')

    types = [
            ('HR',      np.int16),
            ('RAdeg',   np.float64),
            ('DEdeg',   np.float64),
            ('pmRA',    np.float32),
            ('pmDE',    np.float32),
            ('Plx',     np.float32),
            ('n_Plx',   'S1'),
            ('Vmag',    np.float32),
            ('n_Vmag',  'S1'),
            ('u_Vmag',  'S1'),
            ('B-V',     np.float32),
            ('u_B-V',   'S1'),
            ('U-B',     np.float32),
            ('u_U-B',   'S1'),
            ('R-I',     np.float32),
            ('n_R-I',   'S1'),
            ('SpType',  'S20'),
            ('n_SpType','S1'),
            ('RadVel',  np.float32),
            ('n_RadVel','S4'),
            ('l_RotVel','S1'),
            ('RotVel',  np.float32),
            ('u_RotVel','S1'),
            ('Dmag',    np.float32),
            ('Sep',     np.float32),
            ('MultID',  'S4'),
            ('MultCnt', np.int16),
        ]

    tmp = list(zip(*types))
    record = np.dtype({'names': tmp[0], 'formats': tmp[1]})

    data = []
    infile = open(inputfile)
    for row in infile:
        hr = int(row[0:4])
        try:
            ra = (int(row[75:77]) + int(row[77:79])/60. + float(row[79:83])/3600.)*15.
        except:
            ra = np.NaN
        try:
            de = int(row[84:86]) + int(row[86:88])/60. + float(row[88:90])/3600.
        except:
            de = np.NaN
        if row[83]=='-':
            de = -de
        vmag = _str_to_float(row[102:107], np.NaN)
        n_vmag = row[107].strip()
        u_vmag = row[108].strip()
        bv   = _str_to_float(row[109:114], np.NaN)
        u_bv = row[114].strip()
        ub   = _str_to_float(row[115:120], np.NaN)
        u_ub = row[120].strip()
        ri   = _str_to_float(row[121:126], np.NaN)
        u_ri = row[126].strip()
        sptype = row[127:147].strip()
        n_sptype = row[147].strip()
        pmra = _str_to_float(row[148:154], np.NaN)*1000
        pmde = _str_to_float(row[154:160], np.NaN)*1000
        n_plx = row[160].strip()
        plx = _str_to_float(row[161:166], np.NaN)
        radv = _str_to_float(row[166:170], np.NaN)
        if len(row)>=174:
            n_radv = row[170:174].strip()
        else:
            n_radv = ''
        if len(row)>=176:
            l_rotv = row[174:176].strip()
        else:
            l_rotv = ''
        if len(row)>=179:
            rotv = _str_to_float(row[176:179], np.NaN)
        else:
            rotv = np.NaN
        if len(row)>=180:
            u_rotv = row[179].strip()
        else:
            u_rotv = ''
        if len(row)>=180:
            dmag = _str_to_float(row[180:184], np.NaN)
            sep  = _str_to_float(row[184:190], np.NaN)
            multid = row[190:194].strip()
            multcnt = _str_to_int(row[194:196], -1)
        else:
            dmag = np.NaN
            sep  = np.NaN
            multid = ''
            multcnt = -1

        item = (hr, ra, de, pmra, pmde, plx, n_plx, 
                vmag, n_vmag, u_vmag, bv, u_bv, ub, u_ub, ri, u_ri,
                sptype, n_sptype, radv, n_radv, l_rotv, rotv, u_rotv,
                dmag, sep, multid, multcnt)
        data.append(item)
    data = np.array(data, dtype=record)

    pri_hdu = fits.PrimaryHDU()
    tbl_hdu = fits.BinTableHDU(data)
    hdu_lst = fits.HDUList([pri_hdu, tbl_hdu])

    outname = 'BSC.fits'
    if os.path.exists(outname):
        os.remove(outname)
    hdu_lst.writeto(outname)


if __name__=='__main__':
    main()
