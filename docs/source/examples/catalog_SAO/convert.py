#!/usr/bin/env python3
import os
import numpy as np
import astropy.io.fits as fits
from stella.catalog.base import _str_to_float

def main():
    inputfile = os.path.join(os.getenv('ASTRO_DATA'), 'catalog/I/131A/sao.dat')

    types = [
            ('SAO',      np.int32),
            ('RAdeg',    np.float64),
            ('DEdeg',    np.float64),
            ('pmRA',     np.float32),
            ('pmDE',     np.float32),
            ('e_pmRA',   np.float32),
            ('e_pmDE',   np.float32),
            ('Vmag',     np.float32),
            ('Pmag',     np.float32),
            ('SpType',   'S3'),
            ('r_Vmag',   np.int16),
            ('r_Pmag',   np.int16),
            ('r_Num',    np.int16),
            ('r_pm',     np.int16),
            ('r_SpType', np.int16),
            ('RAdeg2000', np.float64),
            ('DEdeg2000', np.float64),
            ('pmRA2000',  np.float32),
            ('pmDE2000',  np.float32),
            ]

    tmp = list(zip(*types))
    record = np.dtype({'names': tmp[0], 'formats': tmp[1]})

    data = []
    infile = open(inputfile)
    for row in infile:
        sao = int(row[0:6])
        radeg = (int(row[7:9]) + int(row[9:11])/60. + float(row[11:17])/3600.)*15.
        dedeg = int(row[42:44]) + int(row[44:46])/60. + float(row[46:51])/3600.
        if row[41]=='-':
            dedeg = -dedeg

        pmra   = _str_to_float(row[17:24], np.NaN)*1000.
        e_pmra = _str_to_float(row[24:26], np.NaN)
        pmde   = _str_to_float(row[51:57], np.NaN)*1000.
        e_pmde = _str_to_float(row[57:59], np.NaN)
        pmag   = _str_to_float(row[76:80], np.NaN)
        vmag   = _str_to_float(row[80:84], np.NaN)
        sptype = row[84:87].strip()
        r_vmag = int(row[87:89])
        r_num  = int(row[89:91])
        r_pmag = int(row[91])
        r_pm   = int(row[92])
        r_sptype = int(row[93])
        ra2000 = (int(row[150:152]) + int(row[152:154])/60. + float(row[154:160])/3600.)*15.
        pmra2000 = _str_to_float(row[160:167], np.NaN)*1000.
        de2000 = int(row[168:170]) + int(row[170:172])/60. + float(row[172:177])/3600.
        pmde2000 = _str_to_float(row[177:183], np.NaN)*1000.
        if row[167]=='-':
            de2000 = -de2000

        item = (sao, radeg, dedeg, pmra, pmde, e_pmra, e_pmde,
                vmag, pmag, sptype, r_vmag, r_pmag, r_num, r_pm, r_sptype,
                ra2000, de2000, pmra2000, pmde2000)

        data.append(item)

    data = np.array(data, dtype=record)

    pri_hdu = fits.PrimaryHDU()
    tbl_hdu = fits.BinTableHDU(data)
    hdu_lst = fits.HDUList([pri_hdu, tbl_hdu])

    outname = 'SAO.fits'
    if os.path.exists(outname):
        os.remove(outname)
    hdu_lst.writeto(outname)



if __name__=='__main__':
    main()
