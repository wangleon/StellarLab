import math

from .errors import ColorIndexError, ParamRangeError, MissingParamError

def color_to_BC(ref, **kwargs):
    '''Get BC using a varietyof calibration relations.

    Notes
    -----
    Available calibration relations:

    * Alonso1995
    * Alonso1999
    * Masana2006
    '''

    if ref == 'Alonso1995':
        bc = _get_dwarf_BC_Alonso1995(**kwargs)
    elif ref == 'Alonso1999':
        bc = _get_giant_BC_Alonso1999(**kwargs)
    elif ref == 'Masana2006':
        bc = _get_dwarf_BC_Masana2006(**kwargs)
    return bc

def _get_dwarf_BC_Alonso1995(**kwargs):
    '''Get BC for dwarfs using the calibration relations given by `Alonso+ 1995
    <http://adsabs.harvard.edu/abs/1995A&A...297..197A>`_.
    
    Notes
    -----
    V-K is in Johnson system

    References
    ----------
    * `Alonso et al., 1995, A&A, 297, 197 <http://adsabs.harvard.edu/abs/1995A&A...297..197A>`_
    '''
    reference = 'Alonso, 1995, A&A, 297, 197'

    teff  = kwargs.pop('teff', None)
    index = kwargs.pop('index', None)
    color = kwargs.pop('color', None)
    FeH   = kwargs.pop('FeH', None)
    band  = kwargs.pop('band', None)
    extrapolation = kwargs.pop('extrapolation',False)

    if FeH == None:
        raise MissingParamError('[Fe/H]', reference)
    if not extrapolation:
        if FeH < -3.0 or FeH > +0.2:
            raise ParamRangeError('[Fe/H]', FeH, reference)
        elif (-0.5 < FeH <= +0.2 and 0.8 < color < 3.0) or \
             (-1.5 < FeH <= -0.5 and 0.9 < color < 2.6) or \
             (-2.5 < FeH <= -1.5 and 1.1 < color < 2.3) or \
             (-3.0 <=FeH <= -2.5 and 1.2 < color < 2.0):
            pass
        else:
            raise ParamRangeError('V-K', color, reference)

    coeff1 = [+2.38619e-4, -1.93659e-4, +6.52621e-5, -7.95862e-6,
              -1.01449e-5, +8.17345e-6, -2.87876e-6, +5.40944e-7]
    coeff2 = [+2.23403e-4, -1.71897e-4, +5.51085e-5, -6.41071e-6,
              -3.71945e-5, +4.99847e-5, -2.41517e-5, +4.10655e-6]

    phi = lambda coeff,color,FeH: \
            coeff[0] + coeff[1]*color + coeff[2]*color**2 + coeff[3]*color**3 +\
            coeff[4]*FeH + coeff[5]*color*FeH + coeff[6]*color**2*FeH +\
            coeff[7]*color**3*FeH

    vk_sun = 1.486
    phi_sun = phi(coeff1,vk_sun,0.0)

    if extrapolation:
        if color <= 1.7:
            bc_v = -2.5*math.log10(phi(coeff1,color,FeH)/phi_sun) - 0.12
            bc_k = -2.5*math.log10(phi(coeff1,color,FeH)/phi_sun) + 1.366
        else:
            bc_v = -2.5*math.log10(phi(coeff2,color,FeH)/phi_sun) - 0.12
            bc_k = -2.5*math.log10(phi(coeff2,color,FeH)/phi_sun) + 1.366
    else:
        if 0.9 < color <= 1.7:
            bc_v = -2.5*math.log10(phi(coeff1,color,FeH)/phi_sun) - 0.12
            bc_k = -2.5*math.log10(phi(coeff1,color,FeH)/phi_sun) + 1.366
        elif 1.7 < color <= 2.9:
            bc_v = -2.5*math.log10(phi(coeff2,color,FeH)/phi_sun) - 0.12
            bc_k = -2.5*math.log10(phi(coeff2,color,FeH)/phi_sun) + 1.366

    if band == 'V':
        return bc_v
    elif band == 'K':
        return bc_k


def _get_giant_BC_Alonso1999(**kwargs):
    '''Get BC for giants using the calibrations relations given by `Alonso+ 1999
    <http://adsabs.harvard.edu/abs/1999A&AS..140..261A>`_.

    References
    ----------
    * `Alonso et al., 1999, A&AS, 140, 261 <http://adsabs.harvard.edu/abs/1999A&AS..140..261A>`_
    '''

    teff = kwargs.pop('teff',None)
    FeH  = kwargs.pop('FeH',0.0)
    extrapolation = kwargs.pop('extrapolation',False)

    logt = math.log10(teff)

    if extrapolation:
        if logt <= 3.66: choose = 17
        else:            choose = 18
    else:
        if   3.50 <= logt <= 3.67 and +0.2 >= FeH > -0.5: choose = 17
        elif 3.56 <= logt <= 3.67 and -0.5 >= FeH > -1.5: choose = 17
        elif 3.58 <= logt <= 3.67 and -1.5 >= FeH > -2.5: choose = 17
        elif 3.61 <= logt <= 3.67 and -2.5 >= FeH > -3.0: choose = 17

        elif 3.65 <= logt <= 3.96 and +0.2 >= FeH > -0.5: choose = 18
        elif 3.65 <= logt <= 3.83 and -0.5 >= FeH > -1.5: choose = 18
        elif 3.65 <= logt <= 3.80 and -1.5 >= FeH > -2.5: choose = 18
        elif 3.65 <= logt <= 3.74 and -2.5 >= FeH > -3.0: choose = 18

        else: raise ValueError

    x = logt - 3.52

    if choose == 17:
        bc = -5.531e-2/x - 0.6177 + 4.420*x - 2.669*x**2 + 0.6943*x*FeH \
             -0.1071*FeH - 8.612e-3*FeH**2
    elif choose == 18:
        bc = -9.930e-2/x + 2.887e-2 + 2.275*x - 4.425*x**2 + 0.3505*x*FeH \
             -5.558e-2*FeH - 5.375e-3*FeH**2
    return bc
    
def _get_dwarf_BC_Masana2006(**kwargs):
    '''Get BC for dwarfs using the calibration relations given by `Masana+ 2006
    <http://adsabs.harvard.edu/abs/2006A&A...450..735M>`_.

    References
    ----------
    * `Masana et al. 2006, A&A, 450, 735 <http://adsabs.harvard.edu/abs/2006A&A...450..735M>`_

    '''
    index = kwargs.pop('index')
    color = kwargs.pop('color')
    FeH   = kwargs.pop('FeH',0.0)
    logg  = kwargs.pop('logg',4.2)
    extrapolation = kwargs.pop('extrapolation',False)

    if index == 'V-K':
        if extrapolation or \
           (-3.0 <  FeH < -1.5 and 1.0  < color < 2.9) or \
           (-1.5 <= FeH < -0.5 and 0.5  < color < 2.9) or \
           (-0.5 <= FeH <  0.0 and 0.4  < color < 3.0) or \
           ( 0.5 <= FeH <  0.5 and 0.35 < color < 2.8):
            if (extrapolation and color < 1.15) or \
               (not extrapolation and
               0.35 < color < 1.15 and 3.25 <= logg <= 4.75):
                bc = 0.1275 + 0.9907*color - 0.0395*color**2 + 0.0693*FeH + \
                     0.0140*FeH**2 + 0.0120*color*FeH - 0.0253*logg

            elif (extrapolation and color >= 1.15) or \
                 (not extrapolation and
                 1.15 <= color < 3.0 and 3.75 <= logg <= 4.75):
                bc = -0.1041 + 1.2600*color - 0.1570*color**2 + 0.1460*FeH + \
                     0.0010*FeH**2 - 0.0631*color*FeH - 0.0079*logg
            else:
                raise ValueError
            return bc
        else:
            raise ValueError
    else:
        raise ValueError
