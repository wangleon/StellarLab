#!/usr/bin/env python3
import os
import astropy.io.fits as fits
import numpy as np
from stella.catalog.base import _str_to_float, _str_to_int

def main():

    astrodata = os.getenv('ASTRO_DATA')

    types = [
            ('EPIC',     np.int32),
            ('Objtype',  'S8'),
            ('Kepflag',  'S3'),
            ('RAdeg',    np.float64),
            ('DEdeg',    np.float64),
            ('pmRA',     np.float32),
            ('pmDE',     np.float32),
            ('e_pmRA',   np.float32),
            ('e_pmDE',   np.float32),
            ('Plx',      np.float32),
            ('e_Plx',    np.float32),
            ('Bmag',     np.float32),
            ('e_Bmag',   np.float32),
            ('Vmag',     np.float32),
            ('e_Vmag',   np.float32),
            ('umag',     np.float32),
            ('e_umag',   np.float32),
            ('gmag',     np.float32),
            ('e_gmag',   np.float32),
            ('rmag',     np.float32),
            ('e_rmag',   np.float32),
            ('imag',     np.float32),
            ('e_imag',   np.float32),
            ('zmag',     np.float32),
            ('e_zmag',   np.float32),
            ('Jmag',     np.float32),
            ('e_Jmag',   np.float32),
            ('Hmag',     np.float32),
            ('e_Hmag',   np.float32),
            ('Kmag',     np.float32),
            ('e_Kmag',   np.float32),
            ('W1mag',    np.float32),
            ('e_W1mag',  np.float32),
            ('W2mag',    np.float32),
            ('e_W2mag',  np.float32),
            ('W3mag',    np.float32),
            ('e_W3mag',  np.float32),
            ('W4mag',    np.float32),
            ('e_W4mag',  np.float32),
            ('kepmag',   np.float32),
            ('Teff',     np.int16),
            ('ue_Teff',  np.int16),
            ('le_Teff',  np.int16),
            ('logg',     np.float32),
            ('ue_logg',  np.float32),
            ('le_logg',  np.float32),
            ('FeH',      np.float32),
            ('ue_FeH',   np.float32),
            ('le_FeH',   np.float32),
            ('Rad',      np.float32),
            ('ue_Rad',   np.float32),
            ('le_Rad',   np.float32),
            ('Mass',     np.float32),
            ('ue_Mass',  np.float32),
            ('le_Mass',  np.float32),
            ('rho',      np.float32),
            ('ue_rho',   np.float32),
            ('le_rho',   np.float32),
            ('Lum',      np.float32),
            ('ue_Lum',   np.float32),
            ('le_Lum',   np.float32),
            ('Dist',     np.float32),
            ('ue_Dist',  np.float32),
            ('le_Dist',  np.float32),
            ('E(B-V)',   np.float32),
            ('ue_E(B-V)',np.float32),
            ('le_E(B-V)',np.float32),
            ]
    
    tmp = list(zip(*types))
    record = np.dtype({'names':tmp[0],'formats':tmp[1]})
    
    fill_item = (0, '', '',
                np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN,
                #ra      dec      pmra    pmdec e_pmra, e_pmdec   plx,   e_plx
                np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN,
                #Bmag             Vmag             umag            gmag
                np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN,
                #rmag             imag              zmag
                np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN,
                # Jmag            Hmag              Kmag
                np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN,
                # w1              w2            w3              w4
                np.NaN,
                # kpmag
                -1, -1, -1, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN,
                # Teff      logg                    FeH
                np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN,
                # Rad                    Mass                   rho
                np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN,
                # Lum                   dist                    EBV
                )

    epic_ranges = {
            1: (201000001, 210000000),
            2: (210000001, 220000000),
            3: (220000001, 230000000),
            4: (230000001, 240000000),
            5: (240000001, 250000000),
            6: (250000001, 251809654),
            }

    for i in range(1,7):
        count = 0
        data = []

        fn = 'epic_%d_19Dec2017.txt'%i
        print(fn)
        filename = os.path.join(astrodata, 'catalog/others/EPIC/%s'%fn)
    
        prev_epic = 0
        infile = open(filename)
        infile.readline()
        _epic = epic_ranges[i][0]
        for row in infile:
            g       = row.split('|')
            epic    = int(g[0])
            objtype = g[6].strip()
            kepflag = g[7].strip()
            ra      = float(g[9])
            dec     = float(g[10])
            pmra    = _str_to_float(g[11], np.NaN)
            e_pmra  = _str_to_float(g[12], np.NaN)
            pmde    = _str_to_float(g[13], np.NaN)
            e_pmde  = _str_to_float(g[14], np.NaN)
            plx     = _str_to_float(g[15], np.NaN)
            e_plx   = _str_to_float(g[16], np.NaN)
            bmag    = _str_to_float(g[17], np.NaN)
            ebmag   = _str_to_float(g[18], np.NaN)
            vmag    = _str_to_float(g[19], np.NaN)
            evmag   = _str_to_float(g[20], np.NaN)
            umag    = _str_to_float(g[21], np.NaN)
            eumag   = _str_to_float(g[22], np.NaN)
            gmag    = _str_to_float(g[23], np.NaN)
            egmag   = _str_to_float(g[24], np.NaN)
            rmag    = _str_to_float(g[25], np.NaN)
            ermag   = _str_to_float(g[26], np.NaN)
            imag    = _str_to_float(g[27], np.NaN)
            eimag   = _str_to_float(g[28], np.NaN)
            zmag    = _str_to_float(g[29], np.NaN)
            ezmag   = _str_to_float(g[30], np.NaN)
            jmag    = _str_to_float(g[31], np.NaN)
            ejmag   = _str_to_float(g[32], np.NaN)
            hmag    = _str_to_float(g[33], np.NaN)
            ehmag   = _str_to_float(g[34], np.NaN)
            kmag    = _str_to_float(g[35], np.NaN)
            ekmag   = _str_to_float(g[36], np.NaN)
            w1mag   = _str_to_float(g[37], np.NaN)
            ew1mag  = _str_to_float(g[38], np.NaN)
            w2mag   = _str_to_float(g[39], np.NaN)
            ew2mag  = _str_to_float(g[40], np.NaN)
            w3mag   = _str_to_float(g[41], np.NaN)
            ew3mag  = _str_to_float(g[42], np.NaN)
            w4mag   = _str_to_float(g[43], np.NaN)
            ew4mag  = _str_to_float(g[44], np.NaN)
            kepmag  = _str_to_float(g[45], np.NaN)
            teff    = _str_to_int(g[46], -1)
            ue_teff = _str_to_int(g[47], -1)
            le_teff = _str_to_int(g[48], -1)
            logg    = _str_to_float(g[49], np.NaN)
            ue_logg = _str_to_float(g[50], np.NaN)
            le_logg = _str_to_float(g[51], np.NaN)
            feh     = _str_to_float(g[52], np.NaN)
            ue_feh  = _str_to_float(g[53], np.NaN)
            le_feh  = _str_to_float(g[54], np.NaN)
            rad     = _str_to_float(g[55], np.NaN)
            ue_rad  = _str_to_float(g[56], np.NaN)
            le_rad  = _str_to_float(g[57], np.NaN)
            mass    = _str_to_float(g[58], np.NaN)
            ue_mass = _str_to_float(g[59], np.NaN)
            le_mass = _str_to_float(g[60], np.NaN)
            rho     = _str_to_float(g[61], np.NaN)
            ue_rho  = _str_to_float(g[62], np.NaN)
            le_rho  = _str_to_float(g[63], np.NaN)
            lum     = _str_to_float(g[64], np.NaN)
            ue_lum  = _str_to_float(g[65], np.NaN)
            le_lum  = _str_to_float(g[66], np.NaN)
            dist    = _str_to_float(g[67], np.NaN)
            ue_dist = _str_to_float(g[68], np.NaN)
            le_dist = _str_to_float(g[69], np.NaN)
            ebv     = _str_to_float(g[70], np.NaN)
            ue_ebv  = _str_to_float(g[71], np.NaN)
            le_ebv  = _str_to_float(g[72], np.NaN)
            
            item = (epic, objtype, kepflag,
                    ra, dec, pmra, pmde, e_pmra, e_pmde, plx, e_plx,
                    bmag, ebmag, vmag, evmag, umag, eumag, gmag, egmag, rmag, ermag,
                    imag, eimag, zmag, ezmag, jmag, ejmag, hmag, ejmag, kmag, ekmag,
                    w1mag, ew1mag, w2mag, ew2mag, w3mag, ew3mag, w4mag, ew4mag, kepmag,
                    teff, ue_teff, le_teff, logg, ue_logg, le_logg,
                    feh,  ue_feh,  le_feh,  rad,  ue_rad,  le_rad,
                    mass, ue_mass, le_mass, rho,  ue_rho,  le_rho,
                    lum,  ue_lum,  le_lum,  dist, ue_dist, le_dist,
                    ebv,  ue_ebv,  le_ebv,
                    )

            if epic > _epic:
                for j in range(epic - _epic):
                    data.append(fill_item)
                    _epic += 1
            data.append(item)

            _epic += 1

            count += 1
            if count % 1000000==0: print(count)

        infile.close()

        print('Data reading finished')
        data = np.array(data, dtype=record)
        print('Begin file writing')

        pri_hdu = fits.PrimaryHDU()
        tbl_hdu = fits.BinTableHDU(data)
        hdu_lst = fits.HDUList([pri_hdu,tbl_hdu])
        outputfile = 'EPIC_%d.fits'%i
        if os.path.exists(outputfile):
            os.remove(outputfile)
        hdu_lst.writeto(outputfile)
    

if __name__=='__main__':
    main()
