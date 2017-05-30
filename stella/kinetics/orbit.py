import math
import numpy as np

def parse_pairwise(arg):
    ''' parse value with error'''
    if (isinstance(arg, list) or isinstance(arg, tuple)) and \
        len(arg)==2:
        return arg
    else:
        raise ValueError

def parse_value_err(arg):
    if isinstance(arg, float) or isinstance(arg, int):
        return arg, None
    elif (isinstance(arg, list) or isinstance(arg, tuple)) and \
        len(arg)==2:
        return arg
    else:
        raise ValueError

def compute_uvw(**kwargs):
    #        ra, dec, rv, parallax, pm_ra, pm_dec,
    #    rv_err=None, parallax_err=None, pm_ra_err=None, pm_dec_err=None,
    #    hand='left'):
    '''compute Galactic UVW for stars
    according to Johnson & Soderblom, 1987, AJ, 93, 864
    Parameters
    -----------
    ra, dec in degree
    pm: proper motion in mas/yr
    rv: radial velocity in km/s
    parallax in mas
    '''
    from astropy.coordinates import SkyCoord
    from ..constant import ALPHA_NGP, DELTA_NGP, L_NCP, AU, tropical_year

    sin = math.sin
    cos = math.cos
    pi  = math.pi

    alpha = ALPHA_NGP/180.*pi
    delta = DELTA_NGP/180.*pi
    theta = L_NCP/180.*pi

    # parse RA and Dec
    if 'eqcoord' in kwargs:
        eqcoord = kwargs.pop('eqcoord')
        if isinstance(eqcoord, SkyCoord):
            icrs = eqcoord.icrs
            ra = icrs.ra.degree
            dec = icrs.dec.degree
        else:
            ra, dec = parse_pairwise(eqcoord)
    elif 'ra' in kwargs and 'dec' in kwargs:
        ra = kwargs.pop('ra')
        dec = kwargs.pop('dec')
    else:
        raise ValueError

    ra  = ra/180.*pi
    dec = dec/180.*pi

    # parse RV
    if 'rv' in kwargs:
        rv, rv_err = parse_value_err(kwargs.pop('rv'))

    # parse distance
    if 'distance' in kwargs:
        d, d_err = parse_value_err(kwargs.pop('distance'))
    elif 'parallax' in kwargs:
        para, para_err = parse_value_err(kwargs.pop('parallax'))
        d = 1000./para
        if para_err is None:
            d_err = None
        else:
            d_err = d*para_err/para
    else:
        raise ValueError

    # parse proper motion
    if 'pm' in kwargs:
        input_pm_ra, input_pm_dec = parse_pairwise(kwargs.pop('pm'))
        pm_ra,  pm_ra_err  = parse_value_err(input_pm_ra)
        pm_dec, pm_dec_err = parse_value_err(input_pm_dec)
        pm_ra  *= 1e-3
        pm_dec *= 1e-3
        if pm_ra_err is not None:
            pm_ra_err  *= 1e-3
        if pm_dec_err is not None:
            pm_dec_err *= 1e-3
    else:
        raise ValueError

    T1 = np.mat([[ cos(theta),  sin(theta),           0],
                 [ sin(theta), -cos(theta),           0],
                 [          0,           0,           1]])
    T2 = np.mat([[-sin(delta),           0,  cos(delta)],
                 [          0,          -1,           0],
                 [ cos(delta),           0, +sin(delta)]])
    T3 = np.mat([[ cos(alpha),  sin(alpha),           0],
                 [ sin(alpha), -cos(alpha),           0],
                 [          0,           0,           1]])

    T = T1*T2*T3

    U_plus = kwargs.pop('U_plus', 'center')
    if U_plus == 'center':
        pass
    elif U_plus == 'anticenter':
        T[0][:] = -T[0][:]
    else:
        raise ValueError


    A1 = np.mat([[  cos(ra),  sin(ra),         0],
                 [  sin(ra), -cos(ra),         0],
                 [        0,        0,        -1]])
    A2 = np.mat([[ cos(dec),        0, -sin(dec)],
                 [        0,       -1,         0],
                 [-sin(dec),        0, -cos(dec)]])
    A = A1*A2

    B = T*A

    k = AU*1e-3/tropical_year/86400  # 1 AU/year in unit of km/s

    x = np.mat([[rv],
                [k*pm_ra*d],
                [k*pm_dec*d]])

    U, V, W = np.array(B*x).flatten()

    if None in [pm_ra_err, pm_dec_err, rv_err, d_err]:
        return UVW(U, V, W, U_plus=U_plus)
    else:
        C = np.mat(np.array(B)**2)

        e11 = rv_err
        e12 = (k*d)**2*(pm_ra_err**2  + (pm_ra*d_err/d)**2)
        e13 = (k*d)**2*(pm_dec_err**2 + (pm_dec*d_err/d)**2)

        e1 = np.mat([[e11],[e12],[e13]])

        e2c = 2.*pm_ra*pm_dec*k**2*d_err**2/d**4

        e = C*e1

        U_err = math.sqrt(e[0,0] + e2c*B[0,1]*B[0,2])
        V_err = math.sqrt(e[1,0] + e2c*B[1,1]*B[1,2])
        W_err = math.sqrt(e[2,0] + e2c*B[2,1]*B[2,2])
        return ((U, U_err), (V, V_err), (W, W_err))

def compute_GalXYZ(**kwargs):
    '''
    comput Galactic position (X, Y, Z)
    '''
    from astropy.coordinates import SkyCoord
    # parse RA and Dec
    if 'eqcoord' in kwargs:
        eqcoord = kwargs.pop('eqcoord')
        if not isinstance(eqcoord, SkyCoord):
            ra, dec = parse_pairwise(eqcoord)
            frame   = kwargs.pop('frame', 'icrs')
            eqcoord = SkyCoord(ra, dec, frame=frame, unit='deg')
        gal = eqcoord.galactic
        l, b = gal.l.degree, gal.b.degree
    elif 'ra' in kwargs and 'dec' in kwargs:
        ra    = kwargs.pop('ra')
        dec   = kwargs.pop('dec')
        frame = kwargs.pop('frame', 'icrs')
        eqcoord = SkyCoord(ra, dec, frame=frame, unit='deg')
        gal = eqcoord.galactic
        l, b = gal.l.degree, gal.b.degree
    elif 'galactic' in kwargs:
        l, b = parse_pairwise(kwargs.pop('galactic'))
    elif 'l' in kwargs and 'b' in kwargs:
        l = kwargs.pop('l')
        b = kwargs.pop('b')
    else:
        raise ValueError

    l = l/180.*math.pi
    b = b/180.*math.pi

    # parse distance
    if 'distance' in kwargs:
        d, d_err = parse_value_err(kwargs.pop('distance'))
    elif 'parallax' in kwargs:
        para, para_err = parse_value_err(kwargs.pop('parallax'))
        d = 1000./para
        if para_err is None:
            d_err = None
        else:
            d_err = d*para_err/para
    else:
        raise ValueError

    R0 = kwargs.pop('R0', 8.5)

    x = R0 - d*math.cos(b)*math.cos(l)
    y = d*math.cos(b)*math.sin(l)
    z = d*math.sin(b)

    return (x, y, z)

def compute_Galorbit(**kwargs):
    potential_lst = kwarargs.pop('potential')

    if 'xyz' in kwargs:
        x, y, z = kwargs.pop('xyz')
        if isinstance(x, list) or isinstance(x, tuple):
            x, y, z = x[0], y[0], z[0]
        elif isinstance(x, float):
            pass
        else:
            raise ValueError

    if 'uvw' in kwargs:
        uvw = kwargs.pop('uvw')
        if isinstance(uvw, UVW):
            u, v, w = UVW.U, UVW.V, UVW.W
            if UVW.U_plus == 'center':
                u = -u
        elif isinstance(uvw, tuple) or isinstance(uvw, list):
            u, _ = parse_value_err(uvw[0])
            v, _ = parse_value_err(uvw[1])
            w, _ = parse_value_err(uvw[2])
        else:
            raise ValueError


