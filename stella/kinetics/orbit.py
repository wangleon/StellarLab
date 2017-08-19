import math
import numpy as np
from astropy.coordinates import SkyCoord

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

def compute_UVW(**kwargs):
    #        ra, dec, rv, parallax, pm_ra, pm_dec,
    #    rv_err=None, parallax_err=None, pm_ra_err=None, pm_dec_err=None,
    #    hand='left'):
    '''Compute Galactic velocity components (*U*, *V*, *W*)

    Parameters
    -----------
    ra : *float*
        Right Ascension in degree at epoch J2000.0
    dec : *float*
        Declination in degree at epoch J2000.0
    eqcoord : *astropy.coordinates.SkyCoord* instance
        Sky coordinate of object. Either (`ra`, `dec`) or `eqcoord` is necessary
    pm : *list* or *tuple*
        Proper motion in mas/yr. Either (`pm_RA`, `pm_Dec`) or ((`pm_RA`,
        `pm_RA_err`), (`pm_Dec`, `pm_Dec_err`))
    rv : *float*, *list* or *tuple*
        Radial velocity in km/s. Either `rv` as a float or (`rv`, `rv_err`)
    parallax : *float*, *list* or *tuple*
        Parallax in mas. Either `parallax` as a float or (`parallax`,
        `parallax_err`)
    U_plus : *string*, [`center`\ \|\ `anticenter`]
        Positive direction (towards Galactic center or anti-center) of *U*
        componenet. Default is `center`

    Returns
    --------
    UVW : *tuple*
        (*U*, *V*, *W*) velocities or ((*U*, *U_err*), (*V*, *V_err*),
        (*W*, *W_err*)) if all the uncertainties to parallax, proper motion and
        radial velocity are given.

    Notes
    -------
    .. |kms| replace:: km s\ :sup:`−1`
    Calculate the Galactic space velocity components (*U*, *V*, *W*) using the
    formula given by `Johnson & Soderblom 1987
    <http://adsabs.harvard.edu/abs/1987AJ.....93..864J>`_.
    The coordinate, parallax, and proper motion are required. The resulting
    velocities are relative to the Sun. The positive direction of *U* is defined
    as towards Galactic center in right-handed system, and towards Galactic
    anticenter in left-handed system. To correct to the local standard of rest
    (LSR), the solar motion (*U*:sub:`LSR`, *V*:sub:`LSR`, *W*:sub:`LSR`) are
    needed, e.g. (9.6 ± 3.9, 14.6 ± 5.0, 9.3 ± 1.0) |kms| (`Reid et al. 2014
    <http://adsabs.harvard.edu/abs/2014ApJ...783..130R>`_).
    

    Examples
    ---------
    Calculate (*U*, *V*, *W*) velocities relative to the sun for HD 9562 (HIP
    7276). The heliocentric radial velocity is −13.3 km/s (`Bensby et al. 2003
    <http://adsabs.harvard.edu/abs/2003A&A...410..527B>`_, Table 2).

    .. code-block:: python

        from stella.catalog.find_catalog import find_HIP
        from stella.kinetics.orbit import compute_UVW

        hip = 7276
        item = find_HIP(hip)
        u, v, w = compute_UVW(ra=item['RAdeg'], dec=item['DEdeg'], parallax=item['Plx'],
                              rv=-13.3, pm=(item['pmRA'], item['pmDE']))
        print('%+6.2f %+6.2f %+6.2f'%(u, v, w))
        # output: -8.86 -26.35 +12.39

    References
    -----------
    * `Bensby et al., 2003, A&A, 410, 527 <http://adsabs.harvard.edu/abs/2003A&A...410..527B>`_
    * `Johnson & Soderblom, 1987, AJ, 93, 864 <http://adsabs.harvard.edu/abs/1987AJ.....93..864J>`_
    * `Reid et al. 2014, ApJ, 783, 130 <http://adsabs.harvard.edu/abs/2014ApJ...783..130R>`_

    '''
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
        return (U, V, W)
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
    Compute Galactic position (*X*, *Y*, *Z*)

    Parameters
    ----------
    ra : *float*
        Right Ascension in degree at epoch J2000.0
    dec : *float*
        Declination in degree at epoch J2000.0
    eqcoord : *astropy.coordinates.SkyCoord* instance
        Sky coordinate of object
    galactic : *list* or *tuple*
        Galactic coordinate (`l`, `b`)
    l : *float*
        Galactic longitude in degree
    b : *float*
        Galactic latitude in degree
    distance : *float*, *list* or *tuple*
        Distance in pc. Either `distance` as a float or (`distance`,
        `distance_err`)
    parallax : *float*, *list* or *tuple*
        Parallax in mas. Either `parallax` as a float or (`parallax`,
        `parallax_err`)
    R0 : *float*
        Solar distance to the Galactic center in kpc

    Returns
    -------


    '''
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

    d *= 1e-3

    x = R0 - d*math.cos(b)*math.cos(l)
    y = d*math.cos(b)*math.sin(l)
    z = d*math.sin(b)

    return (x, y, z)

def compute_Galorbit(**kwargs):
    '''
    Calculate the orbit in the Milky Way

    Parameters
    -----------
    potential : *list*
        List of Galactic potentials
    xyz : *tuple* or *list*
        Galactic positions
    uvw : *tuple* or *list*
        Galactic space velocity
    solar_uvw : *tuple* or *list*
        Solar space velocity
    t : *list*
        List of integration time

    Returns
    --------
    x_lst : *numpy.array*
    y_lst : *numpy.array*
    z_lst : *numpy.array*

    Examples
    ---------
    Calculate the orbit of the Sun

    .. code-block:: python

        solar_uvw = (9.6, 255.2, 9.3) # from Reid et al. 2014
        t_lst = np.arange(0, 0.4, 0.0001) # in Gyr

        x_lst, y_lst, z_lst = orbit.compute_Galorbit(
                                potential = potential_lst,
                                xyz=(R0,0.,0.),
                                uvw=(0.,0.,0.),
                                solar_uvw=solar_uvw,
                                t=t_lst)

    Calculate the orbit of `HD 122563
    <http://simbad.u-strasbg.fr/simbad/sim-id?Ident=HD+122563>`_ (HIP 68594)

    .. code-block:: python

        from stella.catalog  import find_catalog
        hip = 68594
        item = find_catalog.find_HIP2(hip)
        ra, dec = item['RAdeg'], item['DEdeg']
        rv = (-26.58, 0.15) # from SIMBAD
        parallax = (item['Plx'], item['e_Plx'])
        pm = ((item['pmRA'], item['e_pmRA']),(item['pmDE'], item['e_pmDE']))
        uvw = orbit.compute_UVW(ra=ra,dec=dec,rv=rv,parallax=parallax,pm=pm,U_plus='center')
        xyz = orbit.compute_GalXYZ(ra=ra,dec=dec,parallax=parallax,R0=R0)
        x1_lst, y1_lst, z1_lst = orbit.compute_Galorbit(
                                    potential = potential_lst,
                                    xyz=xyz,
                                    uvw=uvw,
                                    solar_uvw=solar_uvw,
                                    t=t_lst)
    '''
    from scipy.integrate import odeint
    from ..constant import pc

    potential_lst = kwargs.pop('potential')

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
        if isinstance(uvw, tuple) or isinstance(uvw, list):
            u, _ = parse_value_err(uvw[0])
            v, _ = parse_value_err(uvw[1])
            w, _ = parse_value_err(uvw[2])
        else:
            raise ValueError

    solar_uvw = kwargs.pop('solar_uvw')
    if isinstance(solar_uvw, tuple) or isinstance(solar_uvw, list):
        solar_u, _ = parse_value_err(solar_uvw[0])
        solar_v, _ = parse_value_err(solar_uvw[1])
        solar_w, _ = parse_value_err(solar_uvw[2])
    else:
        raise ValueError

    target_u = u + solar_u
    target_v = v + solar_v
    target_w = w + solar_w

    t = kwargs.pop('t')
    t_lst = t*1e9*365.2422*86400 # convert Gyr to second

    def derive(var, t, potential_lst):
        x, y, z, vx, vy, vz = var
        acce_lst = np.array([potential.acce_cartesian(x, y, z)
                             for potential in potential_lst])
        ax = acce_lst[:,0].sum()
        ay = acce_lst[:,1].sum()
        az = acce_lst[:,2].sum()
        return [vx/pc, vy/pc, vz/pc, ax, ay, az]

    vx, vy, vz = -target_u, target_v, target_w
    var0 = x, y, z, vx, vy, vz
    sol = odeint(derive, var0, t_lst, args=(potential_lst,))

    x_lst = sol[:,0]
    y_lst = sol[:,1]
    z_lst = sol[:,2]
    r_lst = np.sqrt(x_lst**2 + y_lst**2 + z_lst**2)

    return x_lst, y_lst, z_lst
