#!/usr/bin/env python
import math

from .errors import ColorIndexError, ParamRangeError, MissingParamError

def get_BC(ref, **kwargs):
    ref = ref.strip().lower()

    if ref in ['masana2006']:
        bc = _get_dwarf_bc_Masana2006(**kwargs)
    elif ref in ['alonso1995']:
        bc = _get_dwarf_bc_Alonso1995(**kwargs)
    elif ref in ['alonso1999']:
        bc = _get_giant_bc_Alonso1999(**kwargs)
    return bc

def _get_dwarf_bc_Masana2006(**kwargs):
    """
    get BC based on calibration of Masana et al. 2006, A&A, 450, 735
    (2006A&A...450..735M)
    """
    index = kwargs.pop('index')
    color = kwargs.pop('color')
    feh   = kwargs.pop('feh',0.0)
    logg  = kwargs.pop('logg',4.2)
    extrapolation = kwargs.pop('extrapolation',False)

    if index == 'V-K':
        if extrapolation or \
           (-3.0 < feh < -1.5 and 1.0 < color < 2.9) or \
           (-1.5 <=feh < -0.5 and 0.5 < color < 2.9) or \
           (-0.5 <=feh <  0.0 and 0.4 < color < 3.0) or \
           ( 0.5 <=feh <  0.5 and 0.35< color < 2.8):
            if (extrapolation and color < 1.15) or \
               (not extrapolation and
               0.35 < color < 1.15 and 3.25 <= logg <= 4.75):
                bc = 0.1275 + 0.9907*color - 0.0395*color**2 + 0.0693*feh + \
                     0.0140*feh**2 + 0.0120*color*feh - 0.0253*logg

            elif (extrapolation and color >= 1.15) or \
                 (not extrapolation and
                 1.15 <= color < 3.0 and 3.75 <= logg <= 4.75):
                bc = -0.1041 + 1.2600*color - 0.1570*color**2 + 0.1460*feh + \
                     0.0010*feh**2 - 0.0631*color*feh - 0.0079*logg
            else:
                raise ValueError
            return bc
        else:
            raise ValueError
    else:
        raise ValueError


def _get_dwarf_bc_Alonso1995(**kwargs):
    """
    get BC for dwarfs based on the relation of Alonso et al., 1995, A&A, 297,
    197 (1995A&A...297..197A).
    V-K is in Johnson system
    """
    reference = 'Alonso, 1995, A&A, 297, 197'

    teff  = kwargs.pop('teff', None)
    index = kwargs.pop('index', None)
    color = kwargs.pop('color', None)
    feh   = kwargs.pop('feh', None)
    band  = kwargs.pop('band', None)
    extrapolation = kwargs.pop('extrapolation',False)

    if feh == None:
        raise MissingParamError('[Fe/H]', reference)
    if not extrapolation:
        if feh < -3.0 or feh > +0.2:
            raise ParamRangeError('[Fe/H]', feh, reference)
        elif (-0.5 < feh <= +0.2 and 0.8 < color < 3.0) or \
             (-1.5 < feh <= -0.5 and 0.9 < color < 2.6) or \
             (-2.5 < feh <= -1.5 and 1.1 < color < 2.3) or \
             (-3.0 <=feh <= -2.5 and 1.2 < color < 2.0):
            pass
        else:
            raise ParamRangeError('V-K', color, reference)

    coeff1 = [+2.38619e-4, -1.93659e-4, +6.52621e-5, -7.95862e-6,
              -1.01449e-5, +8.17345e-6, -2.87876e-6, +5.40944e-7]
    coeff2 = [+2.23403e-4, -1.71897e-4, +5.51085e-5, -6.41071e-6,
              -3.71945e-5, +4.99847e-5, -2.41517e-5, +4.10655e-6]

    phi = lambda coeff,color,feh:
            coeff[0] + coeff[1]*color + coeff[2]*color**2 + coeff[3]*color**3 +\
            coeff[4]*feh + coeff[5]*color*feh + coeff[6]*color**2*feh +\
            coeff[7]*color**3*feh

    vk_sun = 1.486
    phi_sun = phi(coeff1,vk_sun,0.0)

    if extrapolation:
        if color <= 1.7:
            bc_v = -2.5*math.log10(phi(coeff1,color,feh)/phi_sun) - 0.12
            bc_k = -2.5*math.log10(phi(coeff1,color,feh)/phi_sun) + 1.366
        else:
            bc_v = -2.5*math.log10(phi(coeff2,color,feh)/phi_sun) - 0.12
            bc_k = -2.5*math.log10(phi(coeff2,color,feh)/phi_sun) + 1.366
    else:
        if 0.9 < color <= 1.7:
            bc_v = -2.5*math.log10(phi(coeff1,color,feh)/phi_sun) - 0.12
            bc_k = -2.5*math.log10(phi(coeff1,color,feh)/phi_sun) + 1.366
        elif 1.7 < color <= 2.9:
            bc_v = -2.5*math.log10(phi(coeff2,color,feh)/phi_sun) - 0.12
            bc_k = -2.5*math.log10(phi(coeff2,color,feh)/phi_sun) + 1.366

    if band == 'V':
        return bc_v
    elif band == 'K':
        return bc_k


def _get_giant_bc_Alonso1999(**kwargs):

    """
    get BC for giants according to Teff and [Fe/H]
    based on calibration of Alonso et al., 1999, A&AS, 140, 261
    (1999A&AS..140..261A)
    """

    teff = kwargs.pop('teff',None)
    feh  = kwargs.pop('feh',0.0)
    extrapolation = kwargs.pop('extrapolation',False)

    logt = math.log10(teff)

    if extrapolation:
        if logt <= 3.66: choose = 17
        else:            choose = 18
    else:
        if   3.50<=logt<=3.67 and +0.2>=feh>-0.5: choose=17
        elif 3.56<=logt<=3.67 and -0.5>=feh>-1.5: choose=17
        elif 3.58<=logt<=3.67 and -1.5>=feh>-2.5: choose=17
        elif 3.61<=logt<=3.67 and -2.5>=feh>-3.0: choose=17

        elif 3.65<=logt<=3.96 and +0.2>=feh>-0.5: choose=18
        elif 3.65<=logt<=3.83 and -0.5>=feh>-1.5: choose=18
        elif 3.65<=logt<=3.80 and -1.5>=feh>-2.5: choose=18
        elif 3.65<=logt<=3.74 and -2.5>=feh>-3.0: choose=18

        else: raise ValueError

    x = logt - 3.52

    if choose == 17:
        bc = -5.531e-2/x - 0.6177 + 4.420*x - 2.669*x**2 + 0.6943*x*feh \
             -0.1071*feh - 8.612e-3*feh**2
    elif choose == 18:
        bc = -9.930e-2/x + 2.887e-2 + 2.275*x - 4.425*x**2 + 0.3505*x*feh \
             -5.558e-2*feh - 5.375e-3*feh**2
    return bc
    
