#!/usr/bin/env python
import math
import numpy as np

class UVW(object):
    '''Galactic UVW velocities
    '''
    def __init__(self, U, V, W, U_plus='center'):
        self.U, self.U_err = parse_value_err(U)
        self.V, self.V_err = parse_value_err(V)
        self.W, self.W_err = parse_value_err(W)
        self.U_plus = U_plus

    def __str__(self):
        return '(U, V, W) = (%f +- %f, %f +- %f, %f +- %f) km/s (U+ = %s)'%(
                self.U, self.U_err, self.V, self.V_err, self.W, self.W_err, self.U_plus)

    def correct_to_LSR(self, ref='Dehnen1998'):
        '''
        Correct solar U, V, W velocity to LSR

        Notes
        ------
        * `Mihalas1968`: Mihalas, D. & Routly, P.M., 1968, *Galactic Astronomy*
            (San Francisco: W.H. Freeman), chap. 5 ("standard solar motion")
        * `Dehnen1998`: Dehnen & Binney, 1998, *MNRAS*, 298, 387
        * `Huang2015`: Huang et al., 2015, *MNRAS*, 449, 162
        '''

        uvw = {
                'Mihalas1968': (10.4, 14.8, 7.3),
                'Dehnen1998': ((10.00, 0.36), ( 5.25, 0.62), (7.17, 0.38)),
                'Huang2015':  (( 7.01, 0.20), (10.13, 0.12), (4.95, 0.09)),
                }
        if ref in uvw:
            u, v, w = uvw[ref]
            if isinstance(u, tuple):
                solar_U, solar_U_err = u
                solar_V, solar_V_err = v
                solar_W, solar_W_err = w
            else:
                solar_U, solar_U_err = u, 0.0
                solar_V, solar_V_err = v, 0.0
                solar_W, solar_W_err = w, 0.0
        else:
            print 'Unknown reference %s'%ref
            raise ValueError

        if self.U_plus == 'center':
            U_LSR = self.U + solar_U
        elif self.U_plus == 'anticenter':
            U_LSR = self.U - solar_U
        else:
            raise ValueError
        V_LSR = self.V + solar_V
        W_LSR = self.W + solar_W

        if None in [self.U_err, self.V_err, self.W_err]:
            return UVW(U_LSR, V_LSR, W_LSR, U_plus = self.U_plus)
        else:
            U_LSR_err = math.sqrt(self.U_err**2 + solar_U_err**2)
            V_LSR_err = math.sqrt(self.V_err**2 + solar_V_err**2)
            W_LSR_err = math.sqrt(self.W_err**2 + solar_W_err**2)
            return UVW((U_LSR, U_LSR_err),
                       (V_LSR, V_LSR_err),
                       (W_LSR, W_LSR_err), U_plus = self.U_plus)


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
        return UVW((U, U_err), (V, V_err), (W, W_err), U_plus=U_plus)

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
        l, b = parse_pairwise(eqcoord)
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

    Rsun = kwargs.pop('Rsun', 8.5e3)

    x = Rsun - d*math.cos(b)*math.cos(l)
    y = d*math.cos(b)*math.sin(l)
    z = d*math.sin(b)

    return (x, y, z)
