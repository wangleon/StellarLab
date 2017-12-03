#!/usr/bin/env python3
import os
import time
import numpy as np
import astropy.io.fits as fits
from stella.parameter.metal import feh_to_z
from stella.evolution.geneva import Geneva

def main():
    n = 500

    geneva = Geneva()

    names = ('logZ', 'mass0', 'logTeff', 'logL', 'age', 'mass')
    formats = (np.float32,np.float32,'(500,)f','(500,)f','(500,)f','(500,)f')
    custom = np.dtype({'names': names, 'formats': formats})

    for logz in np.arange(-3.5, -0.5+1e-6, 0.02):
        t1 = time.time()
        z = 10**logz
        data = []
        for mass0 in np.arange(0.8, 5.0+1e-6, 0.02):
            track = geneva.get_track(mass0=mass0, z=z, n=n)
            record = np.array((logz, mass0, track[0], track[1], track[2], track[3]),
                             dtype=custom)
            data.append(record)
        data = np.array(data, dtype=custom)

        pri_hdu = fits.PrimaryHDU()
        tbl_hdu = fits.BinTableHDU(data)
        hdu_lst = fits.HDUList([pri_hdu, tbl_hdu])
        filename = 'Geneva_grid/logZ%+4.2f.fits'%logz
        if os.path.exists(filename):
            os.remove(filename)
        hdu_lst.writeto(filename)
        t2 = time.time()
        print('%4.2f %10.6f'%(logz, t2-t1))


if __name__=='__main__':
    main()
