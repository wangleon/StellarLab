import math
import numpy as np
import numpy.polynomial as poly

from .errors import ColorIndexError, ParamRangeError, MissingParamError

def get_BC(**kwargs):
    '''Get bolometric correction (BC) using a variety of calibration relations.

    Available calibration relations:

        * `Alonso1995`: returns *BC* in *V* and *K* bands using (*V* − *K*) and
          [Fe/H] for dwarfs.
        * `Alonso1999`: returns *BC* using *T*:sub:`eff` and [Fe/H] for giants.
        * `Masana2006`:
    '''

    ref = kwargs.pop('ref', None)
    if ref == 'Alonso1995':
        bc = _get_dwarf_BC_Alonso1995(**kwargs)
    elif ref == 'Alonso1999':
        bc = _get_giant_BC_Alonso1999(**kwargs)
    elif ref == 'Masana2006':
        bc = _get_dwarf_BC_Masana2006(**kwargs)
    return bc

def _get_BC_Flower1996(teff):
    '''Get *BC* in *V* band according to *T*:sub:`eff` using the relation given
    by `Flower 1996 <http://adsabs.harvard.edu/abs/1996ApJ...469..355F>`_.

    Args:
        teff (integer or float): Effective temperature (*T*:sub:`eff`).
    Returns:
        float: *BC* in *V* band.

    The coefficients given in Table 6 of `Flower 1996
    <http://adsabs.harvard.edu/abs/1996ApJ...469..355F>`_ missed powers of ten.
    `Torres 2010 <http://adsabs.harvard.edu/abs/2010AJ....140.1158T>`_ gave the
    corrent version in Table 1.
    '''
    coeff1 = [-0.118115450538963E+06,
               0.137145973583929E+06,
              -0.636233812100225E+05,
               0.147412923562646E+05,
              -0.170587278406872E+04,
               0.788731721804990E+02,
             ]
    coeff2 = [-0.370510203809015E+05,
               0.385672629965804E+05,
              -0.150651486316025E+05,
               0.261724637119416E+04,
              -0.170623810323864E+03,
             ]
    coeff3 = [-0.190537291496456E+05,
               0.155144866764412E+05,
              -0.421278819301717E+04,
               0.381476328422343E+03,
             ]
    logt = math.log10(teff)

    if logTeff >= 3.9:
        coeff = coeff1
    elif logTeff >= 3.7:
        coeff = coeff2
    else:
        coeff = coeff3

    p = poly.Polynomial(coeff)
    return p(logt)


def _get_dwarf_BC_Alonso1995(**kwargs):
    '''Get *BC* in *V* or *K* for dwarfs using the calibration relations given
    by `Alonso+ 1995 <http://adsabs.harvard.edu/abs/1995A&A...297..197A>`_.
    
    Parameters:
        V_K (float): (*V* − *K*) color
        FeH (float): [Fe/H] ratio
        band (string, optional): Either "V" or "K"

    Returns:
        float or dict: *BC*:sub:`V` or *BC*:sub:`K`, if `band` is given; or
        (*BC*:sub:`V`, *BC*:sub:`K`), if `band` is not given

    Notes:
        The empirical zero points of the Sun are adopted in Johnson system:

        * (*V* −  *K*)\ :sub:`⊙` = 1.486
        * *BC*:sub:`⊙`\ (*V*) = −0.12
        * *BC*:sub:`⊙`\ (*K*) = 1.366

    Examples:

    .. code-block:: python

        from stella.parameter.bc import get_BC
        # find BC in V band
        bc_v, bc_k = get_BC(V_K=1.733, FeH=-0.22, ref='Alonso1995')
        # or
        bc_v, bc_k = get_BC(index='V-K', color=1.733, FeH=-0.22, ref='Alonso1995')

    References:
        * `Alonso et al., 1995, A&A, 297, 197 <http://adsabs.harvard.edu/abs/1995A&A...297..197A>`_
    '''
    reference = 'Alonso, 1995, A&A, 297, 197'

    V_K           = kwargs.pop('V_K', None)
    if V_K is None:
        index = kwargs.pop('index', None)
        if index is not None and index == 'V-K':
            V_K  = kwargs.pop('color', None)

    FeH           = kwargs.pop('FeH', None)
    band          = kwargs.pop('band', None)
    extrapolation = kwargs.pop('extrapolation',False)

    if FeH == None:
        raise MissingParamError('[Fe/H]', reference)
    if not extrapolation:
        if FeH < -3.0 or FeH > +0.2:
            raise ParamRangeError('[Fe/H]', FeH, reference)
        elif (-0.5 < FeH <= +0.2 and 0.8 < V_K < 3.0) or \
             (-1.5 < FeH <= -0.5 and 0.9 < V_K < 2.6) or \
             (-2.5 < FeH <= -1.5 and 1.1 < V_K < 2.3) or \
             (-3.0 <=FeH <= -2.5 and 1.2 < V_K < 2.0):
            pass
        else:
            raise ParamRangeError('V-K', V_K, reference)

    # coefficients coming from equation 9
    coeff1 = np.array([+2.38619e-4, -1.93659e-4, +6.52621e-5, -7.95862e-6,
                       -1.01449e-5, +8.17345e-6, -2.87876e-6, +5.40944e-7])
    # coefficients coming from equation 10
    coeff2 = np.array([+2.23403e-4, -1.71897e-4, +5.51085e-5, -6.41071e-6,
                       -3.71945e-5, +4.99847e-5, -2.41517e-5, +4.10655e-6])

    # coefficients coming from equation 9
    coeff1 = np.array([[+2.38619e-4, -1.93659e-4, +6.52621e-5, -7.95862e-6],
                       [-1.01449e-5, +8.17345e-6, -2.87876e-6, +5.40944e-7]])
    # coefficients coming from equation 10
    coeff2 = np.array([[+2.23403e-4, -1.71897e-4, +5.51085e-5, -6.41071e-6],
                       [-3.71945e-5, +4.99847e-5, -2.41517e-5, +4.10655e-6]])

    phi = lambda coeff: poly.polynomial.polyval2d(FeH, V_K, coeff)

    VK_sun = 1.486
    phi_sun = poly.polynomial.polyval(VK_sun, coeff1[0])

    if extrapolation:
        if V_K <= 1.7:
            bc_v = -2.5*math.log10(phi(coeff1)/phi_sun) - 0.12
            bc_k = -2.5*math.log10(phi(coeff1)/phi_sun) + 1.366
        else:
            bc_v = -2.5*math.log10(phi(coeff2)/phi_sun) - 0.12
            bc_k = -2.5*math.log10(phi(coeff2)/phi_sun) + 1.366
    else:
        if 0.9 < V_K <= 1.7:
            bc_v = -2.5*math.log10(phi(coeff1)/phi_sun) - 0.12
            bc_k = -2.5*math.log10(phi(coeff1)/phi_sun) + 1.366
        elif 1.7 < V_K <= 2.9:
            bc_v = -2.5*math.log10(phi(coeff2)/phi_sun) - 0.12
            bc_k = -2.5*math.log10(phi(coeff2)/phi_sun) + 1.366

    if band is None:
        return (bc_v, bc_k)
    elif band == 'V':
        return bc_v
    elif band == 'K':
        return bc_k
    else:
        return None


def _get_giant_BC_Alonso1999(**kwargs):
    '''Get BC for giants using the calibrations relations given by `Alonso+ 1999
    <http://adsabs.harvard.edu/abs/1999A&AS..140..261A>`_.

    Args:
        Teff (float or int): *T*:sub:`eff` of the star
        FeH (float): [Fe/H] abundance ratio
        extrapolation (bool): use extrapolation of True
    Returns:
        float: *BC*:sub:`V` for the star

    References:
        * `Alonso et al., 1999, A&AS, 140, 261 <http://adsabs.harvard.edu/abs/1999A&AS..140..261A>`_
    '''

    teff          = kwargs.pop('Teff', None)
    FeH           = kwargs.pop('FeH', 0.0)
    extrapolation = kwargs.pop('extrapolation', False)

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
