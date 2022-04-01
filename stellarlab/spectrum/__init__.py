from . import specwcs
from .speclib import get_phoenix_r10kspec, get_phoenix_hires

def wl_vac_to_air(wl_vac, unit='Angstrom' ,ref='Ciddor1996'):
    """
    Convert vacuum wavelength to air wavelength.

    Args:
        wl_vac (float): Vacuum wavelength in unit of Angstrom
        ref (str): Reference.
    Returns:
        float: Air wavelength
    See also:
        * :func:`wl_air_to_vac`
    """

    if ref in ['IAU', 'Edlen1953']:
        a, b1, b2, c1, c2 = 6.4328e-5, 2.94981e-2, 2.5540e-4, 146.0, 41.0
        eq = 1
    elif ref=='Edlen1966':
        a, b1, b2, c1, c2 = 8.34213e-5, 2.406030e-2, 1.5997e-4, 130.0, 38.9
        eq = 1
    elif ref=='Peck1972':
        a, b1, b2, c1, c2 = 0.0, 5.791817e-2, 1.67909e-3, 238.0185, 57.362
        eq = 1
    elif ref=='Birch1994':
        a, b1, b2, c1, c2 = 8.34254e-5, 2.406147e-2, 1.5998e-4, 130.0, 38.9
        eq = 1
    elif ref=='Ciddor1996':
        a, b1, b2, c1, c2 = 0.0, 5.792105e-2, 1.67917e-3, 238.0185, 57.362
        eq = 1
    elif ref=='SDSS':
        eq = 2

    if eq == 1:
        # convert wavelength to wavenumber in micron
        if unit == 'Angstrom':
            wn = 1e4/wl_vac
        elif unit == 'micron':
            wn = 1./wl_vac
        elif unit == 'nm':
            wn = 1e3/wl_vac
        n = 1. + a + b1/(c1-wn**2) + b2/(c2-wn**2)

    elif eq == 2:
        # convert wavelength to wavenumber in Angstrom
        if unit == 'Angstrom':
            wn = 1/wl_vac
        elif unit == 'micron':
            wn = 1e-4/wl_vac
        elif unit == 'nm':
            wn = 1e-1/wl_vac

        n = 1. + 2.735182e-4 + 131.4182*wn**2 + 2.76249e8*wn**4

    return wl_vac/n

def wl_air_to_vac(wl_air, unit='Angstrom' ,ref='Ciddor1996'):
    """
    Convert air wavelength to vacuum wavelength.

    Args:
        wl_air (float): Air wavelength in unit of Angstrom
        ref (str): Reference.
    Returns:
        float: Air wavelength
    See also:
        * :func:`wl_vac_to_air`
    """

    wl = wl_air
    for i in range(3):
        r = wl/wl_vac_to_air(wl, unit=unit, ref=ref)
        wl = wl_air*r
    return wl

def load_stis_x1d(filename):
    """load HST x1d fits
    """
    head = fits.getheader(filename)
    if filename.lower()[-9:]=='_x1d.fits' and \
       head['TELESCOP'] == 'HST' and \
       head['INSTRUME'] == 'STIS':
        data = pf.getdata(filename)
        order_lst = data['SPORDER']
        spec = Spec()
        for o in range(order_lst.size):
            wv = data['WAVELENGTH'][o]
            flux = data['FLUX'][o]
            spd = SpecData(wv=wv, flux=flux)
            spec.add_order(order_lst[o],spd)
        return spec
    else:
        raise ValueError

def load_iue_mxhi(filename):
    '''
    load IUE .MXHI fits
    '''
    head = pf.getheader(filename)
    if filename.lower()[-5:]=='.mxhi' and \
        head['TELESCOP'] == 'IUE':
        data = pf.getdata(filename)
        spec = Spec()
        for i in range(data['ORDER'].size):
            order      = data['ORDER'][i]
            npoints    = data['NPOINTS'][i]
            wavelength = data['WAVELENGTH'][i]
            startpix   = data['STARTPIX'][i]
            deltaw     = data['DELTAW'][i]
            abs_cal    = data['ABS_CAL'][i]
            wv   = np.arange(npoints)*deltaw + wavelength
            flux = abs_cal[startpix-1:startpix-1+npoints]
            spd = SpecData(wv=wv,flux=flux)
            spec.add_order(order,spd)
        return spec
    else:
        raise ValueError

def load_harps_e2ds(filename):
    '''
    load HARPS e2ds fits
    '''
    head = pf.getheader(filename)
    if 'e2ds' in filename.lower() and \
        head['TELESCOP'] == 'ESO-3P6' and\
        head['INSTRUME'] == 'HARPS':
        data = pf.getdata(filename)
        spec = Spec()
        npix  = head['NAXIS1']
        nap   = head['NAXIS2']
        nblue = head['HIERARCH ESO DRS CAL TH ORDER NBLUE']
        nred  = head['HIERARCH ESO DRS CAL TH ORDER NRED']
        ngap  = head['HIERARCH ESO DRS CAL TH ORDER NGAP']
        o1    = head['HIERARCH ESO DRS CAL TH ORDER START']
        ndeg  = head['HIERARCH ESO DRS CAL TH DEG LL']
        for ap in np.arange(nap):
            if ap < nblue:
                order = o1 - ap
            else:
                order = o1 - ap - ngap

            wv = np.zeros(npix)
            for i in np.arange(ndeg+1):
                a = ap*(ndeg+1) + i
                key = 'HIERARCH ESO DRS CAL TH COEFF LL%d'%a
                coeff = head[key]
                wv += coeff*(np.arange(npix)+1)**i
            flux = data[ap]

            spec.add_order(order,SpecData(wv=wv,flux=flux))

        return spec
    else:
        raise ValueError

def load_uves_swca(filename):
    '''
    load UVES SWCA fits
    '''
    head = pf.getheader(filename)
    if head['TELESCOP'] == 'ESO-VLT-U2' and \
       head['INSTRUME'] == 'UVES' and \
       head['NAXIS'] == 2 and \
       head['CTYPE2'] == 'ORDER':
        data = pf.getdata(filename)
        spec = Spec()
        npts   = head['NAXIS1']
        norder = head['NAXIS2']
        for i in range(norder):
            o = i + 1
            wv0 = head['WSTART%d'%o]
            wv1 = head['WEND%d'%o]
            wv = np.linspace(wv0,wv1,npts)
            flux = data[i,:]
            specdata = SpecData(wv=wv,flux=flux)
            spec.add_order(o,specdata)

        return spec

    else:
        raise ValueError

