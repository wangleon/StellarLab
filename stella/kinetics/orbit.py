#!/usr/bin/env python
import math
import numpy as np

class UVW:
    '''Galactic UVW velocities
    uvw = UVW(2.3, 4.5, 5.2, 0.3, 2.1, 1.4, hand='left')
    '''
    def __init__(self, U, V, W,
        U_err=None, V_err=None, W_err=None, hand='left'):
        self.U     = U
        self.V     = V
        self.W     = W
        self.U_err = U_err
        self.V_err = V_err
        self.W_err = W_err
        hand = hand.strip().lower()
        if hand in ['left','right']:
            self.hand = hand
        else:
            raise ValueError

    def __str__(self):
        return '(U, V, W) = (%f +- %f, %f +- %f, %f +- %f) km/s'%(
                self.U, self.U_err, self.V, self.V_err, self.W, self.W_err)

    def correct_to_LSR(self, ref='Dehnen1998'):

        if ref == 'Dehnen1998':
            # solar U, V, W velocity relative to LSR
            # Dehnen & Binney, 1998, MNRAS, 298, 387
            # U in left hand
            solar_U, solar_U_err = -10.0, 0.36
            solar_V, solar_V_err = +5.25, 0.62
            solar_W, solar_W_err = +7.17, 0.38
            solar_hand = 'left'
        elif ref == 'Mihalas1968':
            # Mihalas, D. & Routly, P.M., 1968, Galactic Astronomy (San Francisco: W.H. Freeman), chap. 5
            solar_U, solar_U_err = -10.4, 0.0
            solar_V, solar_V_err = +14.8, 0.0
            solar_W, solar_W_err = +7.3, 0.0
            solar_hand = 'left'

        else:
            print 'Unknown reference %s'%ref
            raise ValueError

        if self.hand == solar_hand:
            U_LSR = self.U + solar_U
        else:
            U_LSR = self.U - solar_U
        V_LSR = self.V + solar_V
        W_LSR = self.W + solar_W

        if None not in [self.U_err, self.V_err, self.W_err]:
            U_LSR_err = math.sqrt(self.U_err**2 + solar_U_err**2)
            V_LSR_err = math.sqrt(self.V_err**2 + solar_V_err**2)
            W_LSR_err = math.sqrt(self.W_err**2 + solar_W_err**2)
        else:
            U_LSR_err, V_LSR_err, W_LSR_err = None, None, None

        return UVW(U_LSR, V_LSR, W_LSR,
                   U_LSR_err, V_LSR_err, W_LSR_err,
                   hand = self.hand)




def compute_uvw(ra, dec, rv, parallax, pm_ra, pm_dec,
        rv_err=None, parallax_err=None, pm_ra_err=None, pm_dec_err=None,
        hand='left'):
    '''compute Galactic UVW for stars
    according to Johnson & Soderblom, 1987, AJ, 93, 864
    ra, dec in degree
    pm: proper motion in mas/yr
    rv: radial velocity in km/s
    parallax in mas
    '''
    from ..constant import ALPHA_NGP, DELTA_NGP, L_NCP
    sin = math.sin
    cos = math.cos
    pi  = math.pi

    pm_ra      *= 1e-3
    pm_dec     *= 1e-3

    alpha = ALPHA_NGP/180.*pi
    delta = DELTA_NGP/180.*pi
    theta = L_NCP/180.*pi


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

    hand = hand.strip().lower()
    if hand == 'left':
        T[0][:] = -T[0][:]
    elif hand == 'right':
        pass
    else:
        raise ValueError

    ra  = ra/180.*pi
    dec = dec/180.*pi

    A1 = np.mat([[  cos(ra),  sin(ra),         0],
                 [  sin(ra), -cos(ra),         0],
                 [        0,        0,        -1]])
    A2 = np.mat([[ cos(dec),        0, -sin(dec)],
                 [        0,       -1,         0],
                 [-sin(dec),        0, -cos(dec)]])
    A = A1*A2

    B = T*A

    k = 4.74057  # 1 AU/year in unit of km/s

    x = np.mat([[rv],
                [k*pm_ra/parallax],
                [k*pm_dec/parallax]])

    U, V, W = np.array(B*x).flatten()

    if None not in [pm_ra_err, pm_dec_err, rv_err, parallax_err]:

        pm_ra_err  *= 1e-3
        pm_dec_err *= 1e-3

        C = np.mat(np.array(B)**2)

        e11 = rv_err
        e12 = (k/parallax)**2*(pm_ra_err**2  +
                                 (pm_ra*parallax_err/parallax)**2  )
        e13 = (k/parallax)**2*( pm_dec_err**2 +
                                 (pm_dec*parallax_err/parallax)**2 )

        e1 = np.mat([[e11],[e12],[e13]])

        e2c = 2.*pm_ra*pm_dec*k**2*parallax_err**2/parallax**4

        e = C*e1

        U_err = math.sqrt(e[0,0] + e2c*B[0,1]*B[0,2])
        V_err = math.sqrt(e[1,0] + e2c*B[1,1]*B[1,2])
        W_err = math.sqrt(e[2,0] + e2c*B[2,1]*B[2,2])
    else:
        U_err, V_err, W_err = None, None, None

    return UVW(U, V, W, U_err, V_err, W_err, hand=hand)
