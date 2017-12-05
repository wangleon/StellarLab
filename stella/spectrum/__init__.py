
def wl_vac_to_air(wl_vac, unit='Angstrom',ref='Ciddor1996'):
    '''
    Convert vacuum wavelength to air wavelength.

    Args:
        wl_vac (float): Vacuum wavelength in unit of Angstrom
        ref (str): Reference.
    Returns:
        float: Air wavelength
    '''

    if ref in ['IAU','Edlen1953']:
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
