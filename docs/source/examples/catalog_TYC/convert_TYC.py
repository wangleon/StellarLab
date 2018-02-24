#!/usr/bin/env python3
import os
import numpy as np
import astropy.io.fits as fits
from stella.catalog.base import _str_to_float, _str_to_int

def main():
    path = os.path.join(os.getenv('ASTRO_DATA'), 'catalog/I/239/')
    
    types = [
            ('TYC',      np.int32),
            ('RAdeg',    np.float64),
            ('DEdeg',    np.float64),
            ('e_RA',     np.float32),
            ('e_DE',     np.float32),
            ('pmRA',     np.float32),
            ('pmDE',     np.float32),
            ('e_pmRA',   np.float32),
            ('e_pmDE',   np.float32),
            ('Plx',      np.float32),
            ('e_Plx',    np.float32),
            ('Nastro',   np.int16),
            ('Vmag',     np.float32),
            ('r_Vmag',   'S1'),
            ('BTmag',    np.float32),
            ('e_BTmag',  np.float32),
            ('VTmag',    np.float32),
            ('e_VTmag',  np.float32),
            ('r_BTmag',  'S1'),
            ('B-V',      np.float32),
            ('e_B-V',    np.float32),
            ('Nphoto',   np.int16),
            ('VTscat',   np.float32),
            ('VTmax',    np.float32),
            ('VTmin',    np.float32),
            ('VarFlag',  'S1'),
            ('MultFlag', 'S1'),
            ]
    tmp = list(zip(*types))
    record = np.dtype({'names':tmp[0],'formats':tmp[1]})

    fill_item = (0, np.NaN, np.NaN, np.NaN, np.NaN, # ra, de, e_ra, e_de
                 np.NaN, np.NaN, np.NaN, np.NaN, # pm
                 np.NaN, np.NaN, -1, # plx, Nastro
                 np.NaN, '', np.NaN, np.NaN, np.NaN, np.NaN, '', # vmag ~ r_BTmag
                 np.NaN, np.NaN, -1, np.NaN, np.NaN, np.NaN, '', ''
                 )
    
    data = []

    fn = 'tyc_main.dat'
    infile = open(os.path.join(path, fn))
    for row in infile:
        if row[0]!='T':
            print('Warning: Not Tycho star')
        g = row[2:14].split()
        tyc1   = int(g[0])
        tyc2   = int(g[1])
        tyc3   = int(g[2])
        tyc = np.int32((tyc1<<18) + (tyc2<<4) + (tyc3<<1))
    
        radeg    = _str_to_float(row[51:63], np.NaN)
        dedeg    = _str_to_float(row[64:76], np.NaN)
        vmag     = _str_to_float(row[41:46], np.NaN)
        r_vmag   = row[49].strip()
        plx      = _str_to_float(row[79:86], np.NaN)
        pm_ra    = _str_to_float(row[87:95], np.NaN)
        pm_de    = _str_to_float(row[96:104], np.NaN)
        e_ra     = _str_to_float(row[105:111], np.NaN)
        e_de     = _str_to_float(row[112:118], np.NaN)
        e_plx    = _str_to_float(row[119:125], np.NaN)
        e_pmra   = _str_to_float(row[126:132], np.NaN)
        e_pmde   = _str_to_float(row[133:139], np.NaN)
        nastro   = _str_to_int(row[200:203], 0)
        btmag    = _str_to_float(row[217:223], np.NaN)
        e_btmag  = _str_to_float(row[224:229], np.NaN)
        vtmag    = _str_to_float(row[230:236], np.NaN)
        e_vtmag  = _str_to_float(row[237:242], np.NaN)
        r_btmag  = row[243].strip()
        bv       = _str_to_float(row[245:251], np.NaN)
        e_bv     = _str_to_float(row[252:257], np.NaN)
        nphoto   = _str_to_int(row[269:272], 0)
        vtscat   = _str_to_float(row[273:278], np.NaN)
        vtmax    = _str_to_float(row[279:284], np.NaN)
        vtmin    = _str_to_float(row[285:290], np.NaN)
        varflag  = row[293].strip()
        multflag = row[295].strip()

        item = (tyc, radeg, dedeg, e_ra, e_de, pm_ra, pm_de, e_pmra, e_pmde,
                plx, e_plx, nastro, 
                vmag, r_vmag, btmag, e_btmag, vtmag, e_vtmag, r_btmag,
                bv, e_bv, nphoto, vtscat, vtmax, vtmin, varflag, multflag)
        data.append(item)
    infile.close()
    
    data = np.array(data, dtype=record)
    data = np.sort(data, order='TYC')
    # add an empty item in the beginning
    data = np.insert(data, 0, np.array(fill_item, dtype=record))
    
    pri_hdu = fits.PrimaryHDU()
    tbl_hdu = fits.BinTableHDU(data)
    hdu_lst = fits.HDUList([pri_hdu,tbl_hdu])
    outputfile = 'TYC.fits'
    if os.path.exists(outputfile):
        os.remove(outputfile)
    hdu_lst.writeto(outputfile)
    
    data = fits.getdata(outputfile)
    tyc_lst = data['TYC']
    diff = np.diff(tyc_lst)
    mask = diff==0
    for i, m in enumerate(mask):
        if m:
            print(data[i]['TYC'], data[i+1]['TYC'])
    
    
if __name__=='__main__':
    main()
