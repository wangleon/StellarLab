#!/usr/bin/env python3
import os
import numpy as np
import astropy.io.fits as fits

def str_to_int(string, exception_value=None):
    try:
        return int(string)
    except:
        return exception_value

def str_to_float(string, exception_value=None):
    try:
        return float(string)
    except:
        return exception_value

def main():
    astrodata = os.getenv('ASTRO_DATA')
    inputfile = os.path.join(astrodata, 'catalog/III/135A/catalog.dat')
    types = [
            ('HD',    np.int32),
            ('RAdeg', np.float64),
            ('DEdeg', np.float64),
            ('q_Ptm', np.int16),
            ('Ptm',   np.float32),
            ('q_Ptg', np.int16),
            ('Ptg',   np.float32),
            ('SpT',   'S3'),
            ('Int',   'S2'),
            ('Rem',   'S1'),
            ]
    tmp = list(zip(*types))
    record = np.dtype({'names': tmp[0], 'formats': tmp[1]})

    data = []
    infile = open(inputfile)
    for row in infile:
        hd    = int(row[0:6])
        rah   = int(row[18:20])
        ram   = int(row[20:23])/10.
        radeg = (rah + ram/60.)*15.
        ded   = int(row[24:26])
        dem   = int(row[26:28])
        dedeg = ded + dem/60.
        if row[23] == '-':
            dedeg = -dedeg
        q_ptm = str_to_int(row[28], -1)
        ptm   = str_to_float(row[29:34], np.NaN)
        q_ptg = str_to_int(row[35], -1)
        ptg   = str_to_float(row[36:41], np.NaN)
        spt   = row[42:45].strip()
        Int   = row[45:47].strip()
        if len(row)<47:
            rem = ''
        else:
            rem = row[47].strip()

        item = np.array((hd, radeg, dedeg, q_ptm, ptm, q_ptg, ptg, spt,
                         Int, rem), dtype=record)
        data.append(item)
    infile.close()

    data = np.array(data, dtype=record)

    pri_hdu = fits.PrimaryHDU()
    tbl_hdu = fits.BinTableHDU(data)
    hdu_lst = fits.HDUList([pri_hdu, tbl_hdu])

    outputfile = 'HD.fits'
    if os.path.exists(outputfile):
        os.remove(outputfile)
    hdu_lst.writeto(outputfile)

if __name__=='__main__':
    main()
