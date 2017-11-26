import numpy.polynomial as poly
from .error import ColorIndexError, ParamRangeError, MissingParamError

def color_to_Teff(index, color, ref, **kwargs):
    '''Convert color to *T*:sub:`eff` using a variety of calibration relations.

    Notes
    -----
    Available calibration ralitions:

    * Flower1996
    * Alonso1996
    * Alonso1999
    * Ramirez2005
    * Masana2006
    * Gonzalez2009
    * Onehag2009
    * Casagrande2010

    References
    ----------
    * `Alonso et al., 1996, A&A, 313, 873 <http://adsabs.harvard.edu/abs/1996A&A...313..873A>`_
    * `Alonso et al., 1999, A&AS, 140, 261 <http://adsabs.harvard.edu/abs/1999A&AS..140..261A>`_
    * `Alonso et al., 2001, A&A, 376, 1039 <http://adsabs.harvard.edu/abs/2001A&A...376.1039A>`_
    * `Casagrande et al., 2010, A&A, 512, 54 <http://adsabs.harvard.edu/abs/2010A&A...512A..54C>`_
    * `Flower, 1996, ApJ, 469, 355 <http://adsabs.harvard.edu/abs/1996ApJ...469..355F>`_
    * `González Hernández & Bonifacio, 2009, A&A, 497, 497 <http://adsabs.harvard.edu/abs/2009A&A...497..497G>`_
    * `Masana et al. 2006, A&A, 450, 735 <http://adsabs.harvard.edu/abs/2006A&A...450..735M>`_
    * `Önehag et al., 2009, A&A, 498, 527 <http://adsabs.harvard.edu/abs/2009A&A...498..527O>`_
    * `Ramírez & Meléndez, 2005, ApJ, 626, 465 <http://adsabs.harvard.edu/abs/2005ApJ...626..465R>`_
    '''
    if ref == 'Flower1996':
        is_supergiant = kwargs.pop('is_supergiant', False)
        Teff = _BV_to_Teff_Flower1996(color, is_supergiant)

    elif ref == 'Alonso1996':
        Teff = _get_dwarf_Teff_Alonso1996(index, color, **kwargs)

    elif ref == 'Alonso1999':
        Teff = _get_giant_Teff_Alonso1999(index, color, **kwargs)

    elif ref == 'Alonso':
        startype = kwargs.pop('startype', None)
        if startype == 'dwarf':
            Teff = _get_dwarf_Teff_Alonso1996(index, color, **kwargs)
        elif startype == 'giant':
            Teff = _get_giant_Teff_Alonso1999(index, color, **kwargs)

    elif ref == 'Ramirez2005':
        startype = kwargs.pop('startype', None)
        if startype == 'dwarf':
            Teff = _get_dwarf_Teff_Ramirez2005(index, color, **kwargs)
        elif startype == 'giant':
            Teff = _get_giant_Teff_Ramirez2005(index, color, **kwargs)

    elif ref == 'Masana2006':
        Teff = _get_dwarf_Teff_Masana2006(index, color, **kwargs)

    elif ref == 'GB2009':
        startype = kwargs.pop('startype', None)
        if startype == 'dwarf':
            Teff = _get_dwarf_Teff_GB2009(index, color, **kwargs)
        elif startype == 'giant':
            Teff = _get_giant_Teff_GB2009(index, color, **kwargs)

    elif ref == 'Onehag2009':
        startype = kwargs.pop('startype', None)
        if startype == 'dwarf':
            Teff = _get_dwarf_Teff_Onehag2009(index, color, **kwargs)
        elif startype == 'giant':
            Teff = _get_giant_teff_Onehag2009(index, color, **kwargs)

    elif ref == 'Casagrande2010':
        Teff = _get_dwarf_Teff_Casagrande2010(index, color, **kwargs)

    return Teff

def _BV_to_Teff_Flower1996(BV, is_supergiant=True):
    '''Convert *B* − *V* to *T*:sub:`eff` using the calibration relations given
    by `Flower 1996 <http://adsabs.harvard.edu/abs/1996ApJ...469..355F>`_.

    Notes
    -----
    The coefficients are given in Table 5 of `Flower 1996
    <http://adsabs.harvard.edu/abs/1996ApJ...469..355F>`_.
    However, the formula and coefficiens were printed incorrectly.
    `Torres 2010 <http://adsabs.harvard.edu/abs/2010AJ....140.1158T>`_ gave the
    correct version in their Table 2.

    References
    ----------
    * `Flower, 1996, ApJ, 469, 355 <http://adsabs.harvard.edu/abs/1996ApJ...469..355F>`_
    * `Torres, 2010, AJ, 140, 1158 <http://adsabs.harvard.edu/abs/2010AJ....140.1158T>`_

    '''
    if is_supergiant:
        coeffs = [4.012559732366214,
                 -1.055043117465989,
                  2.133394538571825,
                 -2.459769794654992,
                  1.349423943497744,
                 -0.283942579112032]
    else:
        coeffs = [3.979145106714099,
                 -0.654992268598245,
                  1.740690042385095,
                 -4.608815154057166,
                  6.792599779944473,
                 -5.396909891322525,
                  2.192970376522490,
                 -0.359495739295671]
    p = poly.Polynomial(coeffs)
    logTeff = p(BV)
    return 10**logTeff

def _get_dwarf_Teff_Alonso1996(index, color, **kwargs):
    '''Convert color and [Fe/H] to *T*:sub:`eff` for dwarfs and subdwarfs using
    the calibration relations given by `Alonso+ 1996
    <http://adsabs.harvard.edu/abs/1996A&A...313..873A>`_.

    References
    ----------
    * `Alonso et al., 1996, A&A, 313, 873 <http://adsabs.harvard.edu/abs/1996A&A...313..873A>`_

    '''

    FeH = kwargs.pop('FeH', 0.0)
    extrapolation = kwargs.pop('extrapolation', False)

    if not extrapolation:
        if FeH < -3.5 or FeH > +0.5:
            raise ParamRangeError('[Fe/H]', FeH, reference)

    if index == 'B-V':
        if not extrapolation:
            if (-0.5 <= FeH <=+0.5 and 0.20 <= color <= 1.5) or \
               (-1.5 <= FeH < -0.5 and 0.30 <= color <= 1.0) or \
               (-2.5 <= FeH < -1.5 and 0.35 <= color <= 0.9) or \
               (-3.5 <= FeH < -2.5 and 0.30 <= color <= 0.8):
                pass
            else:
                raise ParamRangeError(index, color, reference)
        theta = 0.541 + 0.533*color + 0.007*color**2 \
                - 0.019*color*FeH - 0.047*FeH \
                - 0.011*FeH**2

    elif index == 'R-I':
        if not extrapolation:
            if (-0.5 <= FeH <=+0.5 and 0.10 <= color <= 1.00) or \
               (-1.5 <= FeH < -0.5 and 0.20 <= color <= 0.65) or \
               (-2.5 <= FeH < -1.5 and 0.25 <= color <= 0.55) or \
               (-3.5 <= FeH < -2.5 and 0.25 <= color <= 0.50):
                pass
            else:
                raise ParamRangeError(index, color, reference)
        theta = 0.522 + 1.178*color - 0.320*color**2 \
                - 0.087*color*FeH + 0.057*FeH + 0.005*FeH**2

    elif index == 'V-R':
        if not extrapolation:
            if (-0.5 <= FeH <=+0.5 and 0.25 <= color <= 1.40) or \
               (-1.5 <= FeH < -0.5 and 0.30 <= color <= 0.75) or \
               (-2.5 <= FeH < -1.5 and 0.40 <= color <= 0.75) or \
               (-3.5 <= FeH < -2.5 and 0.40 <= color <= 0.70):
                pass
            else:
                raise ParamRangeError(index, color, reference)
        if color <= 0.6:
            theta = 0.474 + 0.755*color + 0.005*color**2 \
                    + 0.003*color*FeH - 0.027*FeH - 0.007*FeH**2
        else:
            theta = 0.524 + 0.724*color - 0.082*color**2 \
                    - 0.166*color*FeH + 0.074*FeH - 0.009*FeH**2

    elif index == 'V-I':
        if not extrapolation:
            if (-0.5 <= FeH <=+0.5 and 0.50 <= color <= 2.50) or \
               (-1.5 <= FeH < -0.5 and 0.60 <= color <= 1.30) or \
               (-2.5 <= FeH < -1.5 and 0.70 <= color <= 1.30) or \
               (-3.5 <= FeH < -2.5 and 0.65 <= color <= 1.20):
                pass
            else:
                raise ParamRangeError(index, color, reference)
        theta = 0.424 + 0.610*color - 0.096*color**2

    elif index == 'V-K':
        if not extrapolation:
            if (-0.5 <= FeH <=+0.5 and 0.4 <= color <= 4.1) or \
               (-1.5 <= FeH < -0.5 and 0.8 <= color <= 3.0) or \
               (-2.5 <= FeH < -1.5 and 1.1 <= color <= 2.4) or \
               (-3.5 <= FeH < -2.5 and 1.1 <= color <= 2.2):
                pass
            else:
                raise ParamRangeError(index, color, reference)
        if color <= 1.6:
            theta = 0.555 + 0.195*color + 0.013*color**2 \
                    - 0.008*color*FeH + 0.009*FeH - 0.002*FeH**2
        else:
            theta = 0.566 + 0.217*color - 0.003*color**2 \
                    - 0.024*color*FeH + 0.037*FeH - 0.002*FeH**2

    elif index == 'b-y':

        # check if c1 exists
        try:
            c1 = kwargs.pop('c1')
        except KeyError:
            raise MissingParamError('c1', reference)

        if not extrapolation:
            if FeH < -3.0 or FeH > +0.5:
                raise ParamRangeError('[Fe/H]', FeH, reference)

            # b-y applicable range according to Figure 11a
            if color < 0.1 or color > 1.0:
                raise ParamRangeError(index, color, reference)

        theta = 0.537 + 0.854*color + 0.196*color**2 \
                -0.198*color*c1 - 0.026*color*FeH \
                -0.014*FeH - 0.009*FeH**2

    elif index == 'beta':
        if not extrapolation:
            if (-0.5 <= FeH <=+0.5 and 2.44 <= color <= 2.74) or \
               (-1.5 <= FeH < -0.5 and 2.50 <= color <= 2.70) or \
               (-2.5 <= FeH < -1.5 and 2.50 <= color <= 2.63) or \
               (-3.5 <= FeH < -2.5 and 2.51 <= color <= 2.62):
                pass
            else:
                raise ParamRangeError(index, color, reference)

        theta = 47.7477 - 34.0506*color + 6.1625*color**2 - 0.1016*color*FeH \
                + 0.3054*FeH + 0.0083*FeH**2

    elif index == 'J-K':
        if not extrapolation:
            if (-0.5 <= FeH <=+0.5 and 0.05 <= color <= 0.85) or \
               (-1.5 <= FeH < -0.5 and 0.15 <= color <= 0.65) or \
               (-2.5 <= FeH < -1.5 and 0.25 <= color <= 0.75) or \
               (-3.5 <= FeH < -2.5 and 0.20 <= color <= 0.60):
                pass
            else:
                raise ParamRangeError(index, color, reference)

        theta = 0.582 + 0.799*color + 0.085*color**2

    elif index == 'J-H':
        if not extrapolation:
            if (-0.5 <= FeH <=+0.5 and 0.00 <= color <= 0.65) or \
               (-1.5 <= FeH < -0.5 and 0.15 <= color <= 0.55) or \
               (-2.5 <= FeH < -1.5 and 0.20 <= color <= 0.60) or \
               (-3.5 <= FeH < -2.5 and 0.15 <= color <= 0.45):
                pass
            else:
                raise ParamRangeError(index, color, reference)
        theta = 0.587 + 0.922*color + 0.218*color**2 + 0.016*color*FeH

    else:
        raise ColorIndexError(index, reference)

    return 5040./theta

def _get_giant_Teff_Alonso1999(index, color, **kwargs):
    '''Convert color and [Fe/H] to *T*:sub:`eff` for giants using the
    calibration relations given by `Alonso+ 1999
    <http://adsabs.harvard.edu/abs/1999A&AS..140..261A>`_.

    References
    ----------
    * `Alonso et al., 1999, A&AS, 140, 261 <http://adsabs.harvard.edu/abs/1999A&AS..140..261A>`_
    * `Alonso et al., 2001, A&A, 376, 1039 <http://adsabs.harvard.edu/abs/2001A&A...376.1039A>`_

    '''

    FeH = kwargs.pop('FeH', 0.0)
    extrapolation = kwargs.pop('extrapolation', False)

    eq = {}
    eq[1] =[0.6388, 0.4065,  -0.1117,   -2.308e-3, -7.783e-2, -1.200e-2]
    eq[2] =[0.8323, 9.374e-2, 1.184e-2,  2.351e-2, -0.1392,   -1.944e-2]
    eq[3] =[0.5716, 0.5404,  -6.126e-2, -4.862e-2, -1.777e-2, -7.969e-3]
    eq[4] =[0.6177, 0.4354,  -4.025e-3,  5.204e-2, -0.1127,   -1.385e-2]
    eq[5] =[0.4972, 0.8841,  -0.1904,   -1.197e-2, -1.025e-2, -5.500e-3]
    #eq[6]=[e-2 (V-I)2-2.693e-2(V-I)3,0.017,125,214
    eq[7] =[0.4974, 1.345,   -0.5008,   -8.134e-2,  3.705e-2, -6.184e-3]
    eq[8] =[0.5558, 0.2105,   1.981e-3, -9.965e-3,  1.325e-2, -2.726e-3]
    eq[9] =[0.3770, 0.3660,  -3.170e-2, -3.074e-3, -2.765e-3, -2.973e-3]
    eq[10]=[0.5977, 1.015,   -1.020e-1, -1.029e-2,  3.006e-2,  1.013e-2]
    eq[11]=[0.5816, 0.9134,  -0.1443,    0.0000,    0.0000,    0.0000  ]
    #eq[12]=[(V-L')*, e-2 (V-L')2-4.651e-3(V-L')3, 0.009, 65, 122
    eq[13]=[0.5859, 0.4846,  -2.457e-2,  0.0000,    0.0000,    0.0000  ]
    eq[14]=[0.5815, 0.7263,   6.856e-2, -6.832e-2, -1.062e-2, -1.079e-2]
    eq[15]=[0.4399, 1.209,   -0.3541,    8.443e-2, -0.1063,   -1.686e-2]
    eq[16]=[0.5883, 0.2008,  -5.931e-3,  5.319e-3, -1.000e-1, -1.542e-2]

    # determine equation
    if index == 'U-V':
        if extrapolation:
            if color <= 1.35: choose = 1
            else:             choose = 2
        else:
            if   0.40<=color<=1.20 and +0.2>=FeH>-0.5: choose = 1
            elif 0.35<=color<=1.20 and -0.5>=FeH>-1.5: choose = 1
            elif 0.40<=color<=1.20 and -1.5>=FeH>-2.5: choose = 1
            elif 0.50<=color<=1.20 and -2.5>=FeH>-3.0: choose = 1

            elif 1.50<=color<=3.50 and +0.2>=FeH>-0.5: choose = 2
            elif 1.50<=color<=3.50 and -0.5>=FeH>-1.5: choose = 2
            elif 1.50<=color<=3.25 and -1.5>=FeH>-2.5: choose = 2

            else: raise ValueError

    elif index == 'B-V':
        if extrapolation:
            if color <= 0.75: choose = 3
            else:             choose = 4
        else:
            if   0.20<=color<=0.80 and +0.2>=FeH>-0.5: choose = 3
            elif 0.35<=color<=0.80 and -0.5>=FeH>-1.5: choose = 3
            elif 0.35<=color<=0.80 and -1.5>=FeH>-2.5: choose = 3
            elif 0.50<=color<=0.80 and -2.5>=FeH>-3.0: choose = 3

            elif 0.70<=color<=1.90 and +0.2>=FeH>-0.5: choose = 4
            elif 0.70<=color<=1.80 and -0.5>=FeH>-1.5: choose = 4
            elif 0.70<=color<=1.35 and -1.5>=FeH>-2.5: choose = 4
            elif 0.70<=color<=1.00 and -2.5>=FeH>-3.0: choose = 4

            else: raise ValueError

    elif index == 'V-R':
        if extrapolation:
            choose = 5
        else:
            if   0.15<=color<=1.70 and +0.2>=FeH>-0.5: choose = 5
            elif 0.45<=color<=1.50 and -0.5>=FeH>-1.5: choose = 5
            elif 0.50<=color<=1.00 and -1.5>=FeH>-2.5: choose = 5
            elif 0.55<=color<=0.85 and -2.5>=FeH>-3.0: choose = 5

            else: raise ValueError

    elif index == 'V-I':
        if extrapolation:
            choose = 6
        else:
            if   0.20<=color<=2.90 and +0.2>=FeH>-0.5: choose = 6
            elif 0.80<=color<=2.00 and -0.5>=FeH>-1.5: choose = 6
            elif 0.85<=color<=2.20 and -1.5>=FeH>-2.5: choose = 6
            elif 1.00<=color<=1.70 and -2.5>=FeH>-3.0: choose = 6

            else: raise ValueError

    elif index == 'R-I':
        if extrapolation:
            choose = 7
        else:
            if   0.15<=color<=1.40 and +0.2>=FeH>-0.5: choose = 7
            elif 0.25<=color<=0.80 and -0.5>=FeH>-1.5: choose = 7
            elif 0.35<=color<=0.70 and -1.5>=FeH>-2.5: choose = 7
            elif 0.40<=color<=0.65 and -2.5>=FeH>-3.0: choose = 7

            else: raise ValueError

    elif index == 'V-K':
        if extrapolation:
            if color <= 2.25: choose = 8
            else:             choose = 9
        else:
            if   0.20<=color<=2.50 and +0.2>=FeH>-0.5: choose = 8
            elif 1.00<=color<=2.50 and -0.5>=FeH>-1.5: choose = 8
            elif 1.20<=color<=2.50 and -1.5>=FeH>-2.5: choose = 8
            elif 1.70<=color<=2.50 and -2.5>=FeH>-3.0: choose = 8

            elif 2.00<=color<=4.90 and +0.2>=FeH>-0.5: choose = 9
            elif 2.00<=color<=4.60 and -0.5>=FeH>-1.5: choose = 9
            elif 2.00<=color<=3.40 and -1.5>=FeH>-2.5: choose = 9
            elif 2.00<=color<=2.80 and -2.5>=FeH>-3.0: choose = 9

            else: raise ValueError

    elif index == 'J-H':
        if extrapolation:
            choose = 10
        else:
            if   0.00<=color<=0.90 and +0.2>=FeH>-0.5: choose = 10
            elif 0.20<=color<=0.80 and -0.5>=FeH>-1.5: choose = 10
            elif 0.30<=color<=0.70 and -1.5>=FeH>-2.5: choose = 10
            elif 0.35<=color<=0.65 and -2.5>=FeH>-3.0: choose = 10

            else: raise ValueError

    elif index == 'J-K':
        if extrapolation:
            choose = 11
        else:
            if   0.00<=color<=1.10 and +0.2>=FeH>-0.5: choose = 11
            elif 0.20<=color<=1.00 and -0.5>=FeH>-1.5: choose = 11
            elif 0.30<=color<=0.90 and -1.5>=FeH>-2.5: choose = 11
            elif 0.40<=color<=0.80 and -2.5>=FeH>-3.0: choose = 11

            else: raise ValueError

    elif index == "V-L'":
        if extrapolation:
            choose = 12
        else:
            if   0.40<=color<=5.00 and +0.2>=FeH>-0.5: choose = 12

            else: raise ValueError

    elif index == '(I-K)J':
        if extrapolation:
            choose = 13
        else:
            if   0.00<=color<=1.90 and +0.2>=FeH>-0.5: choose = 13
            elif 0.50<=color<=1.60 and -0.5>=FeH>-1.5: choose = 13
            elif 0.70<=color<=1.50 and -1.5>=FeH>-2.5: choose = 13
            elif 0.80<=color<=1.20 and -2.5>=FeH>-3.0: choose = 13

            else: raise ValueError

    elif index == 'b-y':
        if extrapolation:
            if color <=0.525: choose = 14
            else:             choose = 15
        else:
            if   0.00<=color<=0.55 and +0.2>=FeH>-0.5: choose = 14
            elif 0.30<=color<=0.55 and -0.5>=FeH>-1.5: choose = 14
            elif 0.35<=color<=0.55 and -1.5>=FeH>-2.5: choose = 14
            elif 0.40<=color<=0.55 and -2.5>=FeH>-3.0: choose = 14

            elif 0.50<=color<=1.00 and +0.2>=FeH>-0.5: choose = 15
            elif 0.50<=color<=0.90 and -0.5>=FeH>-1.5: choose = 15
            elif 0.50<=color<=0.80 and -1.5>=FeH>-2.5: choose = 15
            elif 0.50<=color<=0.70 and -2.5>=FeH>-3.0: choose = 15

            else: raise ValueError

    elif index == 'u-b':
        if extrapolation:
            choose = 16
        else:
            if   1.60<=color<=4.00 and +0.2>=FeH>-0.5: choose = 16
            elif 1.60<=color<=3.70 and -0.5>=FeH>-1.5: choose = 16
            elif 1.60<=color<=3.40 and -1.5>=FeH>-2.5: choose = 16
            elif 1.60<=color<=2.60 and -2.5>=FeH>-3.0: choose = 16

            else: raise ValueError

    else:
        raise ValueError

    if choose == 6:
        theta = 0.5379 + 0.3981*color + 4.432e-2*color**2 - 2.693e-2*color**3
    elif choose == 12:
        theta = 0.5641 + 0.1882*color + 1.890e-2*color**2 - 4.651e-3*color**3
    else:
        a = eq[choose]
        theta =  a[0] + a[1]*color + a[2]*color**2 + a[3]*color*FeH + \
                 a[4]*FeH + a[5]*FeH**2

    return 5040./theta

def _get_dwarf_Teff_Ramirez2005(index, color, **kwargs):
    '''Convert color and [Fe/H] to *T*:sub:`eff` for dwarfs using the
    calibration relations given by `Ramirez+ 2005
    <http://adsabs.harvard.edu/abs/2005ApJ...626..465R>`_.

    References
    ----------
    * `Ramírez & Meléndez, 2005, ApJ, 626, 465 <http://adsabs.harvard.edu/abs/2005ApJ...626..465R>`_

    '''
    reference = 'Ramirez et al., 2005, ApJ, 626, 465'

    try:
        FeH = kwargs.pop('FeH')
    except KeyError:
        raise MissingParamError('[Fe/H]', reference)

    extrapolation = kwargs.pop('extrapolation',False)

    # coeffs in table 2
    coef = {}
    coef['B-V']    = [0.5002,0.6440,-0.0690,-0.0230,-0.0566,-0.0170]
    coef['b-y']    = [0.4129,1.2570,-0.2268,-0.0242,-0.0464,-0.0200]
    coef['Y-V']    = [0.0644,1.7517,-0.5264,-0.0044,-0.0407,-0.0132]
    coef['V-S']    = [0.2417,1.3653,-0.3823,-0.0387,-0.0105,-0.0077]
    coef['B2-V1']  = [0.6019,0.7663,-0.0713,-0.0339,-0.0382,-0.0137]
    coef['B2-G']   = [0.8399,0.4909,-0.0666,-0.0360,-0.0468,-0.0124]
    coef['t']      = [0.7696,0.5927, 0.3439,-0.0437,-0.0143,-0.0088]
    coef['V-RC']   = [0.4333,1.4399,-0.5419,-0.0481,-0.0239,-0.0125]
    coef['V-IC']   = [0.3295,0.9516,-0.2290,-0.0316, 0.0003,-0.0081]
    coef['RC-IC']  = [0.2919,2.1141,-1.0723,-0.0756, 0.0267,-0.0041]
    coef['C42-45'] = [0.5153,0.5963,-0.0572,-0.0573,-0.0221,-0.0018]
    coef['C42-48'] = [0.1601,0.4533,-0.0135,-0.0471, 0.0305,-0.0020]
    coef['BT-VT']  = [0.5619,0.4462,-0.0029, 0.0003,-0.0746,-0.0190]
    coef['V-J2']   = [0.4050,0.4792,-0.0617,-0.0392, 0.0401,-0.0023]
    coef['V-H2']   = [0.4931,0.3056,-0.0241,-0.0396, 0.0678, 0.0020]
    coef['V-K2']   = [0.4942,0.2809,-0.0180,-0.0294, 0.0444,-0.0008]
    coef['VT-K2']  = [0.4886,0.2773,-0.0195,-0.0300, 0.0467,-0.0008]

    if index in coef:
        a = coef[index]
    else:
        raise ColorIndexError(index, reference)

    theta = a[0] + a[1]*color + a[2]*color**2 + a[3]*color*FeH \
            + a[4]*FeH + a[5]*FeH**2

    if index == 'B-V':
        P1 = [-261.548, 684.977, -470.049, 79.8977]
        P2 = [-324.033, 1516.44, -2107.37, 852.150]
        P3 = [30.5985, -46.7882]
        P4 = [139.965, -292.329]
        if extrapolation:
            if   FeH > -0.5: P = P1
            elif FeH > -1.5: P = P2
            elif FeH > -2.5: P = P3
            else:            P = P4
        else:
            if   -0.5 < FeH < +0.5 and 0.310 <= color <= 1.507: P = P1
            elif -1.5 < FeH <=-0.5 and 0.307 <= color <= 1.202: P = P2
            elif -2.5 < FeH <=-1.5 and 0.335 <= color <= 1.030: P = P3
            elif -4.0 < FeH <=-2.5 and 0.343 <= color <= 0.976: P = P4
            else:
                raise ParamRangeError(index, color, reference)
    elif index == 'b-y':
        P1 = [-1237.11, 6591.29, -11061.3, 5852.18]
        P2 = [-2617.66, 22607.4, -68325.4, 86072.5, -38602.2]
        P3 = [103.927, -312.419, 225.430]
        P4 = [-294.106, 648.320]
        if extrapolation:
            if   FeH > -0.5: P = P1
            elif FeH > -1.5: P = P2
            elif FeH > -2.5: P = P3
            else:            P = P4
        else:
            if   -0.5 < FeH < +0.5 and 0.248 <= color <= 0.824: P = P1
            elif -1.5 < FeH <=-0.5 and 0.234 <= color <= 0.692: P = P2
            elif -2.5 < FeH <=-1.5 and 0.290 <= color <= 0.672: P = P3
            elif -4.0 < FeH <=-2.5 and 0.270 <= color <= 0.479: P = P4
            else:
                raise ParamRangeError(index, color, reference)
    elif index == 'Y-V':
        P1 = [-10407.1, 42733.6, -27378.8, -96466.3, 162033.0, -70956.4]
        P2 = [11.6451]
        P3 = [-507.732, 1943.73, -1727.66]
        P4 = [-310.166, 496.709]
        if extrapolation:
            if   FeH > -0.5: P = P1
            elif FeH > -1.5: P = P2
            elif FeH > -2.5: P = P3
            else:            P = P4
        else:
            if   -0.5 < FeH < +0.5 and 0.420 <= color <= 0.940: P = P1
            elif -1.5 < FeH <=-0.5 and 0.452 <= color <= 0.660: P = P2
            elif -2.5 < FeH <=-1.5 and 0.455 <= color <= 0.720: P = P3
            elif -4.0 < FeH <=-2.5 and 0.446 <= color <= 0.643: P = P4
            else:
                raise ParamRangeError(index, color, reference)
    elif index == 'V-S':
        P1 = [-1436.48, 5566.00, -6780.53, 2613.40]
        P2 = [-728.818, 2256.18, -1704.54]
        P3 = [101.031, 114.354, -447.778]
        P4 = [596.461, -1130.92]
        if extrapolation:
            if   FeH > -0.5: P = P1
            elif FeH > -1.5: P = P2
            elif FeH > -2.5: P = P3
            else:            P = P4
        else:
            if   -0.5 < FeH < +0.5 and 0.370 <= color <= 1.130: P = P1
            elif -1.5 < FeH <=-0.5 and 0.410 <= color <= 0.690: P = P2
            elif -2.5 < FeH <=-1.5 and 0.441 <= color <= 0.810: P = P3
            elif -4.0 < FeH <=-2.5 and 0.438 <= color <= 0.584: P = P4
            else:
                raise ParamRangeError(index, color, reference)
    elif index == 'B2-V1':
        P1 = [-439.817, 2637.06, -4762.80, 2606.79]
        P2 = [-257.527, 2078.96, -4919.04, 3685.65, -500.348]
        P3 = [-28.5544, 228.735, -295.958]
        P4 = [64.2911, -365.124]
        if extrapolation:
            if   FeH > -0.5: P = P1
            elif FeH > -1.5: P = P2
            elif FeH > -2.5: P = P3
            else:            P = P4
        else:
            if   -0.5 < FeH < +0.5 and 0.119 <= color <= 0.936: P = P1
            elif -1.5 < FeH <=-0.5 and 0.132 <= color <= 0.593: P = P2
            elif -2.5 < FeH <=-1.5 and 0.178 <= color <= 0.621: P = P3
            elif -4.0 < FeH <=-2.5 and 0.185 <= color <= 0.435: P = P4
            else:
                raise ParamRangeError(index, color, reference)
    elif index == 'B2-G':
        P1 = [-6.29800, 160.976, -386.520, 250.628]
        P2 = [21.3254, -56.4562, -651.533, 720.639]
        P3 = [11.5114, -34.5752, -265.563]
        P4 = [-156.547, -313.408, 4886.53]
        if extrapolation:
            if   FeH > -0.5: P = P1
            elif FeH > -1.5: P = P2
            elif FeH > -2.5: P = P3
            else:            P = P4
        else:
            if   -0.5 < FeH < +0.5 and -0.271 <= color <= 1.110: P = P1
            elif -1.5 < FeH <=-0.5 and -0.262 <= color <= 0.502: P = P2
            elif -2.5 < FeH <=-1.5 and -0.200 <= color <= 0.544: P = P3
            elif -4.0 < FeH <=-2.5 and -0.179 <= color <= 0.150: P = P4
            else:
                raise ParamRangeError(index, color, reference)
    elif index == 't':
        P1 = [-16.3530, 273.725, -1383.02, 2274.81]
        P2 = [35.2419, -185.953]
        P3 = [11.9635]
        P4 = [-39.1918]
        if extrapolation:
            if   FeH > -0.5: P = P1
            elif FeH > -1.5: P = P2
            elif FeH > -2.5: P = P3
            else:            P = P4
        else:
            if   -0.5 < FeH < +0.5 and -0.119 <= color <= 0.450: P = P1
            elif -1.5 < FeH <=-0.5 and -0.066 <= color <= 0.373: P = P2
            elif -2.5 < FeH <=-1.5 and -0.006 <= color <= 0.333: P = P3
            elif -4.0 < FeH <=-2.5 and  0.020 <= color <= 0.295: P = P4
            else:
                raise ParamRangeError(index, color, reference)
    elif index == 'V-RC':
        P1 = [-2666.55, 27264.5, -103923.0, 174663.0, -104940.0, -23249.4, 32644.9]
        P2 = [4.20153]
        P3 = [123.940, -342.217]
        P4 = [8.55498]
        if extrapolation:
            if   FeH > -0.5: P = P1
            elif FeH > -1.5: P = P2
            elif FeH > -2.5: P = P3
            else:            P = P4
        else:
            if   -0.5 < FeH < +0.5 and 0.204 <= color <= 0.880: P = P1
            elif -1.5 < FeH <=-0.5 and 0.284 <= color <= 0.546: P = P2
            elif -2.5 < FeH <=-1.5 and 0.264 <= color <= 0.532: P = P3
            elif -4.0 < FeH <=-2.5 and 0.240 <= color <= 0.336: P = P4
            else:
                raise ParamRangeError(index, color, reference)
    elif index == 'V-IC':
        P1 = [-2757.79, 9961.33, -10546.6, -1746.05, 10512.3, -6653.57, 1301.21]
        P2 = [-22.9008, 40.2078]
        P3 = [-667.732, 1709.88, -1069.62]
        if extrapolation:
            if   FeH > -0.5: P = P1
            elif FeH > -1.5: P = P2
            else:            P = P3
        else:
            if   -0.5 < FeH < +0.5 and 0.491 <= color <= 1.721: P = P1
            elif -1.5 < FeH <=-0.5 and 0.597 <= color <= 1.052: P = P2
            elif -2.5 < FeH <=-1.5 and 0.547 <= color <= 1.026: P = P3
            else:
                raise ParamRangeError(index, color, reference)
    elif index == 'RC-IC':
        P1 = [-3326.97, 26263.8, -75355.8, 94246.5, -43334.8]
        P2 = [12.4740]
        P3 = [-5837.31, 41439.2, -94729.8, 69584.8]
        P4 = [32.1826]
        if extrapolation:
            if   FeH > -0.5: P = P1
            elif FeH > -1.5: P = P2
            elif FeH > -2.5: P = P3
            else:            P = P4
        else:
            if   -0.5 < FeH < +0.5 and 0.242 <= color <= 0.838: P = P1
            elif -1.5 < FeH <=-0.5 and 0.300 <= color <= 0.718: P = P2
            elif -2.5 < FeH <=-1.5 and 0.283 <= color <= 0.551: P = P3
            elif -4.0 < FeH <=-2.5 and 0.290 <= color <= 0.364: P = P4
            else:
                raise ParamRangeError(index, color, reference)
    elif index == 'C42-45':
        P1 = [1533.40, -5546.94, 6324.29, -2254.52]
        P2 = [808.065, -2725.54, 2806.13, -902.995]
        if extrapolation:
            if   FeH > -0.5: P = P1
            else:            P = P2
        else:
            if   -0.5 < FeH < +0.5 and 0.461 <= color <= 1.428: P = P1
            elif -1.5 < FeH <=-0.5 and 0.480 <= color <= 0.812: P = P2
            else:
                raise ParamRangeError(index, color, reference)
    elif index == 'C42-48':
        P1 = [658.568, -283.310, -709.877, 575.693, -114.834]
        P2 = [176.678, -204.699, 53.2421]
        P3 = [1069.18, -678.907]
        if extrapolation:
            if   FeH > -0.5: P = P1
            elif FeH > -1.5: P = P2
            else:            P = P3
        else:
            if   -0.5 < FeH < +0.5 and 1.286 <= color <= 2.711: P = P1
            elif -1.5 < FeH <=-0.5 and 1.465 <= color <= 1.957: P = P2
            elif -2.5 < FeH <=-1.5 and 1.399 <= color <= 1.509: P = P3
            else:
                raise ParamRangeError(index, color, reference)
    elif index == 'BT-VT':
        P1 = [1199.21, -5470.57, 8367.46, -5119.55, 1078.09]
        P2 = [-64.1045, 140.575, -59.4233]
        P3 = [-6030.19, 29153.4, -25882.7, -64112.9, 126115.0, -59817.9]
        P4 = [-3255.07, 16259.9, -20315.3]
        if extrapolation:
            if   FeH > -0.5: P = P1
            elif FeH > -1.5: P = P2
            elif FeH > -2.5: P = P3
            else:            P = P4
        else:
            if   -0.5 < FeH < +0.5 and 0.344 <= color <= 1.715: P = P1
            elif -1.5 < FeH <=-0.5 and 0.391 <= color <= 1.556: P = P2
            elif -2.5 < FeH <=-1.5 and 0.380 <= color <= 0.922: P = P3
            elif -4.0 < FeH <=-2.5 and 0.367 <= color <= 0.504: P = P4
            else:
                raise ParamRangeError(index, color, reference)
    elif index == 'V-J2':
        P1 = [422.406, -910.603, 621.335, -132.566]
        P2 = [-466.616, 658.349, -220.454]
        P3 = [-862.072, 1236.84, -423.729]
        P4 = [-1046.10, 1652.06, -597.340]
        if extrapolation:
            if   FeH > -0.5: P = P1
            elif FeH > -1.5: P = P2
            elif FeH > -2.5: P = P3
            else:            P = P4
        else:
            if   -0.5 < FeH < +0.5 and 0.815 <= color <= 2.608: P = P1
            elif -1.5 < FeH <=-0.5 and 0.860 <= color <= 2.087: P = P2
            elif -2.5 < FeH <=-1.5 and 0.927 <= color <= 1.983: P = P3
            elif -4.0 < FeH <=-2.5 and 0.891 <= color <= 1.932: P = P4
            else:
                raise ParamRangeError(index, color, reference)
    elif index == 'V-H2':
        P1 = [-53.5574, 36.0990, 15.6878, -8.84468]
        P2 = [1.60629]
        P3 = [506.559, -1277.52, 939.519, -208.621]
        P4 = [-471.588, 643.972, -199.639]
        if extrapolation:
            if   FeH > -0.5: P = P1
            elif FeH > -1.5: P = P2
            elif FeH > -2.5: P = P3
            else:            P = P4
        else:
            if   -0.5 < FeH < +0.5 and 0.839 <= color <= 3.215: P = P1
            elif -1.5 < FeH <=-0.5 and 1.032 <= color <= 2.532: P = P2
            elif -2.5 < FeH <=-1.5 and 1.070 <= color <= 2.535: P = P3
            elif -4.0 < FeH <=-2.5 and 1.093 <= color <= 2.388: P = P4
            else:
                raise ParamRangeError(index, color, reference)
    elif index == 'V-K2':
        P1 = [-1425.36, 3218.36, -2566.54, 859.644, -102.554]
        P2 = [2.35133]
        P3 = [-1849.46, 4577.00, -4284.02, 1770.38, -268.589]
        P4 = [215.721, -796.519, 714.423, -175.678]
        if extrapolation:
            if   FeH > -0.5: P = P1
            elif FeH > -1.5: P = P2
            elif FeH > -2.5: P = P3
            else:            P = P4
        else:
            if   -0.5 < FeH < +0.5 and 0.896 <= color <= 3.360: P = P1
            elif -1.5 < FeH <=-0.5 and 1.060 <= color <= 2.665: P = P2
            elif -2.5 < FeH <=-1.5 and 1.101 <= color <= 2.670: P = P3
            elif -4.0 < FeH <=-2.5 and 1.126 <= color <= 2.596: P = P4
            else:
                raise ParamRangeError(index, color, reference)
    elif index == 'VT-K2':
        P1 = [-1581.85, 3273.10, -2395.38, 736.352, -80.8177]
        P2 = [68.1279, -130.968, 52.8391]
        P3 = [-2384.82, 4196.14, -2557.04, 595.365, -31.9955]
        P4 = [-628.682, 423.682]
        if extrapolation:
            if   FeH > -0.5: P = P1
            elif FeH > -1.5: P = P2
            elif FeH > -2.5: P = P3
            else:            P = P4
        else:
            if   -0.5 < FeH < +0.5 and 0.942 <= color <= 3.284: P = P1
            elif -1.5 < FeH <=-0.5 and 1.078 <= color <= 2.561: P = P2
            elif -2.5 < FeH <=-1.5 and 1.237 <= color <= 2.406: P = P3
            elif -4.0 < FeH <=-2.5 and 1.170 <= color <= 1.668: P = P4
            else:
                raise ParamRangeError(index, color, reference)
    else:
        raise ColorIndexError(index, reference)

    Teff = 5040./theta

    poly = P[-1]
    for i in range(len(P)-1):
        # P=[x,x,x,x], i=0, 1, 2, 3
        # res = ((P[3]*x+P[2])*x+P[1])*x+P[0]
        poly = poly*color + P[-(i+2)]
    Teff += poly

    return Teff

def _get_giant_Teff_Ramirez2005(index, color, **kwargs):
    '''Convert color and [Fe/H] to *T*:sub:`eff` for giants using the
    calibration relations given by `Ramirez+ 2005
    <http://adsabs.harvard.edu/abs/2005ApJ...626..465R>`_.

    Parameters
    ----------
    
    Notes
    -----
    Avaiable colors and

    References
    ----------
    * `Ramírez & Meléndez, 2005, ApJ, 626, 465 <http://adsabs.harvard.edu/abs/2005ApJ...626..465R>`_
    
    '''
    reference = 'Ramirez et al., 2005, ApJ, 626, 465'

    try:
        FeH = kwargs.pop('FeH')
    except KeyError:
        raise MissingParamError('[Fe/H]', reference)

    extrapolation = kwargs.pop('extrapolation',False)

    # coeffs in table 3
    coef = {}
    coef['B-V']    = [0.5737,0.4882,-0.0149, 0.0563,-0.1160,-0.0114]
    coef['b-y']    = [0.5515,0.9085,-0.1494, 0.0616,-0.0668,-0.0083]
    coef['Y-V']    = [0.3672,1.0467,-0.1995, 0.0650,-0.0913,-0.0133]
    coef['V-S']    = [0.3481,1.1188,-0.2068, 0.0299,-0.0481,-0.0083]
    coef['B2-V1']  = [0.6553,0.6278,-0.0629, 0.0627,-0.0816,-0.0084]
    coef['B2-G']   = [0.8492,0.4344,-0.0365, 0.0466,-0.0696,-0.0107]
    coef['t']      = [0.7460,0.8151,-0.1943, 0.0855,-0.0421,-0.0034]
    coef['V-RC']   = [0.3849,1.6205,-0.6395, 0.1060,-0.0875,-0.0089]
    coef['V-IC']   = [0.3575,0.9069,-0.2025, 0.0395,-0.0551,-0.0061]
    coef['RC-IC']  = [0.4351,1.6549,-0.7215,-0.0610, 0.0332,-0.0023]
    coef['C42-45'] = [0.4783,0.7748,-0.1361,-0.0712,-0.0117, 0.0071]
    coef['C42-48'] = [0.0023,0.6401,-0.0632,-0.0023,-0.0706,-0.0070]
    coef['BT-VT']  = [0.5726,0.4461,-0.0324, 0.0518,-0.1170,-0.0094]
    coef['V-J2']   = [0.2943,0.5604,-0.0677, 0.0179,-0.0532,-0.0088]
    coef['V-H2']   = [0.4354,0.3405,-0.0263,-0.0012,-0.0049,-0.0027]
    coef['V-K2']   = [0.4405,0.3272,-0.0252,-0.0016,-0.0053,-0.0040]
    coef['VT-K2']  = [0.4813,0.2871,-0.0203,-0.0045, 0.0062,-0.0019]

    if index in coef.keys():
        a = coef[index]
    else:
        raise ColorIndexError(index, reference)

    theta = a[0] + a[1]*color + a[2]*color**2 + a[3]*color*FeH \
            + a[4]*FeH + a[5]*FeH**2

    if index == 'B-V':
        P = [112.116, -372.622, 67.1254, 395.333, -203.471]
        P = [-12.9762]
        P = [606.032, -1248.79, 627.453]
        P = [-9.26209]
        if extrapolation:
            if   FeH > -0.5: P = P1
            elif FeH > -1.5: P = P2
            elif FeH > -2.5: P = P3
            else:            P = P4
        else:
            if   -0.5 < FeH < +0.5 and 0.144 <= color <= 1.668: P = P1
            elif -1.5 < FeH <=-0.5 and 0.664 <= color <= 1.558: P = P2
            elif -2.5 < FeH <=-1.5 and 0.605 <= color <= 1.352: P = P3
            elif -4.0 < FeH <=-2.5 and 0.680 <= color <= 1.110: P = P4
            else:
                raise ParamRangeError(index, color, reference)
    elif index == 'b-y':
        P = [-124.159, 553.827, -490.703]
        P = [888.088, -2879.23, 2097.89]
        P = [1867.63, -6657.49, 5784.81]
        P = [348.237, -659.093]
        if extrapolation:
            if   FeH > -0.5: P = P1
            elif FeH > -1.5: P = P2
            elif FeH > -2.5: P = P3
            else:            P = P4
        else:
            if   -0.5 < FeH < +0.5 and 0.053 <= color <= 1.077: P = P1
            elif -1.5 < FeH <=-0.5 and 0.309 <= color <= 0.893: P = P2
            elif -2.5 < FeH <=-1.5 and 0.388 <= color <= 0.702: P = P3
            elif -4.0 < FeH <=-2.5 and 0.404 <= color <= 0.683: P = P4
            else:
                raise ParamRangeError(index, color, reference)
    elif index == 'Y-V':
        P1 = [-308.851, 1241.57, -1524.60, 593.157]
        P2 = [-36.6533, 383.901, -458.085]
        P3 = [3038.83, -8668.15, 6067.04]
        P4 = [2685.88, -7433.07, 4991.81]
        if extrapolation:
            if   FeH > -0.5: P = P1
            elif FeH > -1.5: P = P2
            elif FeH > -2.5: P = P3
            else:            P = P4
        else:
            if   -0.5 < FeH < +0.5 and 0.230 <= color <= 1.290: P = P1
            elif -1.5 < FeH <=-0.5 and 0.558 <= color <= 0.940: P = P2
            elif -2.5 < FeH <=-1.5 and 0.544 <= color <= 0.817: P = P3
            elif -4.0 < FeH <=-2.5 and 0.510 <= color <= 0.830: P = P4
            else:
                raise ParamRangeError(index, color, reference)
    elif index == 'V-S':
        P1 = [-1605.54, 9118.16, -17672.6, 14184.1, -4023.76]
        P2 = [187.841, -270.092]
        P3 = [10.1750]
        P4 = [-14.2019]
        if extrapolation:
            if   FeH > -0.5: P = P1
            elif FeH > -1.5: P = P2
            elif FeH > -2.5: P = P3
            else:            P = P4
        else:
            if   -0.5 < FeH < +0.5 and 0.261 <= color <= 1.230: P = P1
            elif -1.5 < FeH <=-0.5 and 0.508 <= color <= 0.992: P = P2
            elif -2.5 < FeH <=-1.5 and 0.529 <= color <= 0.990: P = P3
            elif -4.0 < FeH <=-2.5 and 0.573 <= color <= 0.790: P = P4
            else:
                raise ParamRangeError(index, color, reference)
    elif index == 'B2-V1':
        P1 = [-15.0383, 50.8876, -32.3978]
        P2 = [80.1344, -147.055]
        P3 = [323.889, -1031.06, 795.024]
        P4 = [1403.86, -4866.09, 4029.75]
        if extrapolation:
            if   FeH > -0.5: P = P1
            elif FeH > -1.5: P = P2
            elif FeH > -2.5: P = P3
            else:            P = P4
        else:
            if   -0.5 < FeH < +0.5 and -0.079 <= color <= 1.321: P = P1
            elif -1.5 < FeH <=-0.5 and  0.385 <= color <= 1.021: P = P2
            elif -2.5 < FeH <=-1.5 and  0.307 <= color <= 0.958: P = P3
            elif -4.0 < FeH <=-2.5 and  0.407 <= color <= 0.648: P = P4
            else:
                raise ParamRangeError(index, color, reference)
    elif index == 'B2-G':
        P1 = [-0.52642, 10.4471, -7.53155]
        P2 = [26.1904, -89.2171]
        P3 = [9.87980]
        P4 = [232.248, -1452.43, 1848.07]
        if extrapolation:
            if   FeH > -0.5: P = P1
            elif FeH > -1.5: P = P2
            elif FeH > -2.5: P = P3
            else:            P = P4
        else:
            if   -0.5 < FeH < +0.5 and -0.543 <= color <= 1.230: P = P1
            elif -1.5 < FeH <=-0.5 and  0.155 <= color <= 0.966: P = P2
            elif -2.5 < FeH <=-1.5 and  0.132 <= color <= 0.991: P = P3
            elif -4.0 < FeH <=-2.5 and  0.104 <= color <= 0.437: P = P4
            else:
                raise ParamRangeError(index, color, reference)
    elif index == 't':
        P = [-46.1506, -60.1641, 643.522, -599.555]
        P = [27.8739, -84.1166]
        P = [67.1191, -139.127]
        P = [122.254, -394.604]
        if extrapolation:
            if   FeH > -0.5: P = P1
            elif FeH > -1.5: P = P2
            elif FeH > -2.5: P = P3
            else:            P = P4
        else:
            if   -0.5 < FeH < +0.5 and 0.072 <= color <= 0.970: P = P1
            elif -1.5 < FeH <=-0.5 and 0.064 <= color <= 0.766: P = P2
            elif -2.5 < FeH <=-1.5 and 0.166 <= color <= 0.619: P = P3
            elif -4.0 < FeH <=-2.5 and 0.215 <= color <= 0.511: P = P4
            else:
                raise ParamRangeError(index, color, reference)
    elif index == 'V-RC':
        P1 = [-8.51797, 15.6675]
        P2 = [-10.7764]
        P3 = [61.9821, -78.7382]
        P4 = [27.9886, -100.149]
        if extrapolation:
            if   FeH > -0.5: P = P1
            elif FeH > -1.5: P = P2
            elif FeH > -2.5: P = P3
            else:            P = P4
        else:
            if   -0.5 < FeH < +0.5 and 0.299 <= color <= 1.106: P = P1
            elif -1.5 < FeH <=-0.5 and 0.387 <= color <= 0.752: P = P2
            elif -2.5 < FeH <=-1.5 and 0.429 <= color <= 0.598: P = P3
            elif -4.0 < FeH <=-2.5 and 0.394 <= color <= 0.550: P = P4
            else:
                raise ParamRangeError(index, color, reference)
    elif index == 'V-IC':
        P1 = [0.42933]
        P2 = [-0.14180]
        P3 = [9.31011]
        P4 = [-23.0514]
        if extrapolation:
            if   FeH > -0.5: P = P1
            elif FeH > -1.5: P = P2
            elif FeH > -2.5: P = P3
            else:            P = P4
        else:
            if   -0.5 < FeH < +0.5 and 0.573 <= color <= 2.000: P = P1
            elif -1.5 < FeH <=-0.5 and 0.795 <= color <= 1.524: P = P2
            elif -2.5 < FeH <=-1.5 and 0.870 <= color <= 1.303: P = P3
            elif -4.0 < FeH <=-2.5 and 0.812 <= color <= 1.095: P = P4
            else:
                raise ParamRangeError(index, color, reference)
    elif index == 'RC-IC':
        P1 = [61.3557, -116.711]
        P2 = [-16.8645]
        P3 = [32.0870]
        P4 = [-15.6318]
        if extrapolation:
            if   FeH > -0.5: P = P1
            elif FeH > -1.5: P = P2
            elif FeH > -2.5: P = P3
            else:            P = P4
        else:
            if   -0.5 < FeH < +0.5 and 0.413 <= color <= 0.793: P = P1
            elif -1.5 < FeH <=-0.5 and 0.383 <= color <= 0.771: P = P2
            elif -2.5 < FeH <=-1.5 and 0.434 <= color <= 0.725: P = P3
            elif -4.0 < FeH <=-2.5 and 0.364 <= color <= 0.545: P = P4
            else:
                raise ParamRangeError(index, color, reference)
    elif index == 'C42-45':
        P = [-68.3798, 109.259, -34.4503]
        P = [-0.62507]
        P = [-40.0150, 35.6803]
        P = [-314.177, 636.443]
        if extrapolation:
            if   FeH > -0.5: P = P1
            elif FeH > -1.5: P = P2
            elif FeH > -2.5: P = P3
            else:            P = P4
        else:
            if   -0.5 < FeH < +0.5 and 0.409 <= color <= 1.369: P = P1
            elif -1.5 < FeH <=-0.5 and 0.430 <= color <= 1.270: P = P2
            elif -2.5 < FeH <=-1.5 and 0.441 <= color <= 0.894: P = P3
            elif -4.0 < FeH <=-2.5 and 0.490 <= color <= 0.640: P = P4
            else:
                raise ParamRangeError(index, color, reference)
    elif index == 'C42-48':
        P1 = [1006.40, -549.012, -649.212, 534.912, -100.038]
        P2 = [-6.92065]
        P3 = [-113.222, 57.3030]
        P4 = [566.914, -329.631]
        if extrapolation:
            if   FeH > -0.5: P = P1
            elif FeH > -1.5: P = P2
            elif FeH > -2.5: P = P3
            else:            P = P4
        else:
            if   -0.5 < FeH < +0.5 and 1.531 <= color <= 2.767: P = P1
            elif -1.5 < FeH <=-0.5 and 1.400 <= color <= 2.647: P = P2
            elif -2.5 < FeH <=-1.5 and 1.466 <= color <= 2.260: P = P3
            elif -4.0 < FeH <=-2.5 and 1.571 <= color <= 1.799: P = P4
            else:
                raise ParamRangeError(index, color, reference)
    elif index == 'BT-VT':
        P1 = [346.881, -1690.16, 2035.65, -797.248, 70.7799]
        P2 = [196.416, -372.164, 126.196]
        P3 = [938.789, -1919.98, 929.779]
        P4 = [1112.46, -2717.81, 1577.18]
        if extrapolation:
            if   FeH > -0.5: P = P1
            elif FeH > -1.5: P = P2
            elif FeH > -2.5: P = P3
            else:            P = P4
        else:
            if   -0.5 < FeH < +0.5 and 0.123 <= color <= 1.953: P = P1
            elif -1.5 < FeH <=-0.5 and 0.424 <= color <= 1.644: P = P2
            elif -2.5 < FeH <=-1.5 and 0.534 <= color <= 1.356: P = P3
            elif -4.0 < FeH <=-2.5 and 0.465 <= color <= 1.026: P = P4
            else:
                raise ParamRangeError(index, color, reference)
    elif index == 'V-J2':
        P = [-122.595, 76.4847]
        P = [-10.3848]
        P = [4.18695, 13.8937]
        P = [-67.7716, 28.9202]
        if extrapolation:
            if   FeH > -0.5: P = P1
            elif FeH > -1.5: P = P2
            elif FeH > -2.5: P = P3
            else:            P = P4
        else:
            if   -0.5 < FeH < +0.5 and 1.259 <= color <= 2.400: P = P1
            elif -1.5 < FeH <=-0.5 and 1.030 <= color <= 3.418: P = P2
            elif -2.5 < FeH <=-1.5 and 1.033 <= color <= 2.679: P = P3
            elif -4.0 < FeH <=-2.5 and 0.977 <= color <= 2.048: P = P4
            else:
                raise ParamRangeError(index, color, reference)
    elif index == 'V-H2':
        P1 = [-377.022, 334.733, -69.8093]
        P2 = [71.7949, -55.5383, 9.61821]
        P3 = [-27.4190, 20.7082]
        P4 = [-46.2946, 20.1061]
        if extrapolation:
            if   FeH > -0.5: P = P1
            elif FeH > -1.5: P = P2
            elif FeH > -2.5: P = P3
            else:            P = P4
        else:
            if   -0.5 < FeH < +0.5 and 1.194 <= color <= 3.059: P = P1
            elif -1.5 < FeH <=-0.5 and 1.293 <= color <= 4.263: P = P2
            elif -2.5 < FeH <=-1.5 and 1.273 <= color <= 3.416: P = P3
            elif -4.0 < FeH <=-2.5 and 1.232 <= color <= 2.625: P = P4
            else:
                raise ParamRangeError(index, color, reference)
    elif index == 'V-K2':
        P1 = [-72.6664, 36.5361]
        P2 = [86.0358, -65.4928, 10.8901]
        P3 = [-6.96153, 14.3298]
        P4 = [-943.925, 1497.64, -795.867, 138.965]
        if extrapolation:
            if   FeH > -0.5: P = P1
            elif FeH > -1.5: P = P2
            elif FeH > -2.5: P = P3
            else:            P = P4
        else:
            if   -0.5 < FeH < +0.5 and 1.244 <= color <= 3.286: P = P1
            elif -1.5 < FeH <=-0.5 and 1.366 <= color <= 4.474: P = P2
            elif -2.5 < FeH <=-1.5 and 1.334 <= color <= 3.549: P = P3
            elif -4.0 < FeH <=-2.5 and 1.258 <= color <= 2.768: P = P4
            else:
                raise ParamRangeError(index, color, reference)
    elif index == 'VT-K2':
        P1 = [-37.2128, 31.2900, -6.72743]
        P2 = [-193.512, 166.183, -33.2781]
        P3 = [-2.02136]
        P4 = [8.06982]
        if extrapolation:
            if   FeH > -0.5: P = P1
            elif FeH > -1.5: P = P2
            elif FeH > -2.5: P = P3
            else:            P = P4
        else:
            if   -0.5 < FeH < +0.5 and 1.107 <= color <= 3.944: P = P1
            elif -1.5 < FeH <=-0.5 and 1.403 <= color <= 3.157: P = P2
            elif -2.5 < FeH <=-1.5 and 1.339 <= color <= 3.750: P = P3
            elif -4.0 < FeH <=-2.5 and 1.668 <= color <= 2.722: P = P4
            else:
                raise ParamRangeError(index, color, reference)
    else:
        raise ColorIndexError(index, reference)

    Teff = 5040./theta

    poly = P[-1]
    for i in range(len(P)-1):
        # P=[x,x,x,x], i=0, 1, 2, 3
        # res = ((P[3]*x+P[2])*x+P[1])*x+P[0]
        poly = poly*color + P[-(i+2)]
    Teff += poly

    return Teff

def _get_dwarf_Teff_Masana2006(index, color, **kwargs):
    '''Convert color, [Fe/H] and log\ *g* to *T*:sub:`eff` for dwarfs using the
    calibration relations given by `Masana+ 2006
    <http://adsabs.harvard.edu/abs/2006A&A...450..735M>`_.

    References
    ----------
    * `Masana et al. 2006, A&A, 450, 735 <http://adsabs.harvard.edu/abs/2006A&A...450..735M>`_

    '''
    reference = 'Masana et al., 2006, A&A, 450, 735'

    if index != 'V-K':
        raise ColorIndexError(index, reference)

    try:
        FeH  = kwargs.pop('FeH')
    except KeyError:
        raise MissingParamError('[Fe/H]', reference)

    try:
        logg = kwargs.pop('logg')
    except KeyError:
        raise MissingParamError('logg', reference)

    extrapolation = kwargs.pop('extrapolation',False)

    if extrapolation or \
       (-3.0 <= FeH < -1.5 and 1.0  <= color <= 2.9) or \
       (-1.5 <= FeH < -0.5 and 0.5  <= color <= 2.9) or \
       (-0.5 <= FeH <  0.0 and 0.4  <= color <= 3.0) or \
       ( 0.5 <= FeH <= 0.5 and 0.35 <= color <= 2.8):
        if (extrapolation and color < 1.15) or \
            (not extrapolation and 0.35<color<1.15 and 3.25<=logg<=4.75):
            theta = 0.5961 + 0.1567*color + 0.0309*color**2 \
                    + 0.009*FeH + 0.0022*FeH**2 \
                    +0.0021*color*FeH - 0.0067*logg
        elif (extrapolation and color>=1.15) or \
            (not extrapolation and 1.15<=color<3.0 and 3.75<=logg<=4.75):
            theta = 0.5135 + 0.2687*color - 0.0174*color**2 \
                    + 0.0298*FeH - 0.0009*FeH**2 \
                    - 0.0184*color*FeH - 0.0028*logg
        else:
            raise ParamRangeError('V-K', color, reference)
        return 5040./theta
    else:
        raise ParamRangeError('V-K', color, reference)

def _get_dwarf_Teff_GB2009(index, color, **kwargs):
    '''Convert color to *T*:sub:`eff` for dwarfs using the calibration relations
    given by `González Hernández & Bonifacio, 2009
    <http://adsabs.harvard.edu/abs/2009A&A...497..497G>`_.

    References
    ----------
    * `González Hernández & Bonifacio, 2009, A&A, 497, 497 <http://adsabs.harvard.edu/abs/2009A&A...497..497G>`_
    '''

    reference = 'Gonzalez et al., 2009, A&A, 497, 497'

    try:
        FeH = kwargs.pop('FeH')
    except KeyError:
        raise MissingParamError('[Fe/H]', reference)

    extrapolation = kwargs.pop('extrapolation',False)

    FeH_range = {}; color_range = {}
    FeH_range['B-V'] = [-3.5,0.5]; color_range['B-V'] = [0.2,1.3]
    FeH_range['V-R'] = [-3.1,0.3]; color_range['V-R'] = [0.2,0.8]
    FeH_range['V-I'] = [-3.1,0.3]; color_range['V-I'] = [0.5,1.4]
    FeH_range['V-J'] = [-3.5,0.5]; color_range['V-J'] = [0.5,2.3]
    FeH_range['V-H'] = [-3.5,0.5]; color_range['V-H'] = [0.6,2.8]
    FeH_range['V-Ks']= [-3.5,0.5]; color_range['V-Ks']= [0.7,3.0]
    FeH_range['J-Ks']= [-3.5,0.5]; color_range['J-Ks']= [0.1,0.8]

    coef = {}
    coef['B-V'] = [0.5725, 0.4722,  0.0086, -0.0628, -0.0038, -0.0051]
    coef['V-R'] = [0.4451, 1.4561, -0.6893, -0.0944,  0.0161, -0.0038]
    coef['V-I'] = [0.4025, 0.8324, -0.2041, -0.0555,  0.0410, -0.0003]
    coef['V-J'] = [0.4997, 0.3504, -0.0230, -0.0295,  0.0468,  0.0037]
    coef['V-H'] = [0.5341, 0.2517, -0.0100, -0.0236,  0.0523,  0.0044]
    coef['V-Ks']= [0.5201, 0.2511, -0.0118, -0.0186,  0.0408,  0.0033]
    coef['J-Ks']= [0.6524, 0.5813,  0.1225, -0.0646,  0.0370,  0.0016]

    if extrapolation:
        b = coef[index]
    else:
        if FeH_range[index][0] <= FeH <= FeH_range[index][1] and \
           color_range[index][0] <= color <= color_range[index][1]:
            b = coef[index]
        else:
            raise ParamRangeError(index, color, reference)

    theta = b[0] + b[1]*color + b[2]*color**2 + b[3]*color*FeH \
            + b[4]*FeH + b[5]*FeH**2
    return 5040./theta

def _get_giant_Teff_GB2009(index, color, **kwargs):
    '''Convert color to *T*:sub:`eff` for giants using the calibration relations
    given by `González Hernández & Bonifacio, 2009
    <http://adsabs.harvard.edu/abs/2009A&A...497..497G>`_.

    References
    ----------
    * `González Hernández & Bonifacio, 2009, A&A, 497, 497 <http://adsabs.harvard.edu/abs/2009A&A...497..497G>`_
    '''

    reference = 'Gonzalez et al., 2009, A&A, 497, 497'

    try:
        FeH = kwargs.pop('FeH')
    except KeyError:
        raise MissingParamError('[Fe/H]', reference)

    extrapolation = kwargs.pop('extrapolation',False)

    FeH_range = {}; color_range = {}
    FeH_range['B-V'] = [-4.0,0.2]; color_range['B-V'] = [0.3,1.4]
    FeH_range['V-R'] = [-4.0,0.1]; color_range['V-R'] = [0.3,0.7]
    FeH_range['V-J'] = [-4.0,0.2]; color_range['V-J'] = [1.0,2.4]
    FeH_range['V-H'] = [-4.0,0.2]; color_range['V-H'] = [0.8,3.1]
    FeH_range['V-Ks']= [-4.0,0.2]; color_range['V-Ks']= [1.1,3.4]
    FeH_range['J-Ks']= [-4.0,0.2]; color_range['J-Ks']= [0.1,0.9]

    coef = {}
    coef['B-V'] = [0.4967, 0.7260, -0.1563,  0.0255, -0.0585, -0.0061]
    coef['V-R'] = [0.4530, 1.4347, -0.5883, -0.0156, -0.0096, -0.0039]
    coef['V-J'] = [0.4629, 0.4124, -0.0417, -0.0012,  0.0094,  0.0013]
    coef['V-H'] = [0.5321, 0.2649, -0.0146, -0.0069,  0.0211,  0.0009]
    coef['V-Ks']= [0.5293, 0.2489, -0.0119, -0.0042,  0.0135,  0.0010]
    coef['J-Ks']= [0.6517, 0.6312,  0.0168, -0.0381,  0.0256,  0.0013]

    if extrapolation:
        b = coef[index]
    else:
        if FeH_range[index][0] <= FeH <= FeH_range[index][1] and \
           color_range[index][0] <= color <= color_range[index][1]:
            b = coef[index]
        else:
            raise ParamRangeError(index, color, reference)

    theta = b[0] + b[1]*color + b[2]*color**2 + b[3]*color*FeH \
            + b[4]*FeH + b[5]*FeH**2
    return 5040./theta

def _get_dwarf_Teff_Onehag2009(index, color, **kwargs):
    '''Convert color and [Fe/H] to *T*:sub:`eff` for dwarfs using the
    calibration relations given by `Önehag+ 2009
    <http://adsabs.harvard.edu/abs/2009A&A...498..527O>`_.

    Notes
    -----
    Only used for stars with log\ *g* >= 4.0

    References
    ----------
    * `Önehag et al., 2009, A&A, 498, 527 <http://adsabs.harvard.edu/abs/2009A&A...498..527O>`_
    '''
    reference = 'Onehag et al., 2009, A&A, 498, 527'

    try:
        FeH = kwargs.pop('FeH')
    except KeyError:
        raise MissingParamError('[Fe/H]', reference)

    extrapolation = kwargs.pop('extrapolation',False)

    if index == 'b-y':

        try:
            c1 = kwargs.pop('c1')
        except KeyError:
            raise MissingParamError('c1', reference)

        if not extrapolation:
            
            if (-3.00 <= Feh <= 0.50 and 0.20 <= color <= 0.70 and \
                 0.10 <= c1 <= 0.55):
                pass
            else:
                raise ParamRangeError(index, color, reference)

        theta = 0.415 + 1.313*color - 0.477*color**2 + 0.188*color*c1 \
                - 0.037*color*FeH - 0.003*FeH - 0.006*FeH**2

    elif index == 'Hbeta':
        if not extrapolation:

            if ((2.44 <= color <= 2.74 and -0.5 <= FeH <= 0.5) or \
                (2.50 <= color <= 2.70 and -1.5 < FeH <= -0.5) or \
                (2.50 <= color <= 2.63 and -2.5 < FeH <= -1.5) or \
                (2.51 <= color <= 2.62 and -3.5 < FeH <= -2.5)):
                pass
            else:
                raise ParamRangeError(index, color, reference)

        theta = 28.60 - 19.79*color + 3.504*color**2 \
                + 0.422*color*FeH - 1.068*FeH + 0.002*FeH**2
    else:
        raise ColorIndexError(index, reference)

    return 5040./theta

def _get_giant_Teff_Onehag2009(index, color, **kwargs):
    '''Convert color and [Fe/H] to *T*:sub:`eff` for giants using the
    calibration relations given by `Önehag+ 2009
    <http://adsabs.harvard.edu/abs/2009A&A...498..527O>`_.

    Notes
    -----
    Only used for 1.5 <= log\ *g* <= 3.5

    References
    ----------
    * `Önehag et al., 2009, A&A, 498, 527 <http://adsabs.harvard.edu/abs/2009A&A...498..527O>`_
    '''
    reference = 'Onehag et al., 2009, A&A, 498, 527'

    try:
        FeH = kwargs.pop('FeH')
    except KeyError:
        raise MissingParamError('[Fe/H]', reference)

    extrapolation = kwargs.pop('extrapolation',False)

    if index == 'b-y':
        coef = {}
        coef[1] = [0.6732,0.0859, 1.1455,-1.080e-2,-0.132e-2,-0.082e-2]
        coef[2] = [0.1983,2.0931,-0.9978, 4.709e-2,-2.66e-2, -0.104e-2]
        coef[3] = [0.4522,1.1745,-0.3093,-0.1693,   2.165e-2,-1.679e-2]
        if extrapolation:
            if color <= 0.424: choose = 1
            elif FeH <= -0.5:  choose = 2
            else:              choose = 3
        else:
            if   0.15  <= color <= 0.424:                         choose = 1
            elif 0.424 <= color <= 0.712 and -5.0 <= FeH <= -0.5: choose = 2
            elif 0.428 <= color <= 0.794 and -0.5 <  FeH <=  0.5: choose = 3
            else:
                raise ParamRangeError(index, color, reference)

        a = coef[choose]
        theta = a[0] + a[1]*color + a[2]*color**2 \
                + a[3]*color*FeH + a[4]*FeH + a[5]*FeH**2
    else:
        raise ColorIndexError(index, reference)

    return 5040./theta

def _get_dwarf_Teff_Casagrande2010(index, color, **kwargs):
    '''Convert color and [Fe/H] to *T*:sub:`eff` for dwarfs using the
    calibration relations given by `Casagrande+ 2010
    <http://adsabs.harvard.edu/abs/2010A&A...512A..54C>`_.

    References
    ----------
    * `Casagrande et al., 2010, A&A, 512, 54 <http://adsabs.harvard.edu/abs/2010A&A...512A..54C>`_
    '''
    reference = 'Casagrande, 2010, A&A, 512, A54'

    try:
        FeH  = kwargs.pop('FeH')
    except KeyError:
        raise MissingParamError('[Fe/H]', reference)

    FeH_range = {}; color_range = {}
    FeH_range['B-V']   = [-5.0,0.4]; color_range['B-V']   = [0.18,1.29]
    FeH_range['V-RC']  = [-5.0,0.3]; color_range['V-RC']  = [0.24,0.80]
    FeH_range['RC-IC'] = [-5.0,0.3]; color_range['RC-IC'] = [0.23,0.68]
    FeH_range['V-IC']  = [-5.0,0.3]; color_range['V-IC']  = [0.46,1.47]
    FeH_range['V-J']   = [-5.0,0.4]; color_range['V-J']   = [0.61,2.44]
    FeH_range['V-H']   = [-5.0,0.4]; color_range['V-H']   = [0.67,3.01]
    FeH_range['V-Ks']  = [-5.0,0.4]; color_range['V-Ks']  = [0.78,3.15]
    FeH_range['J-Ks']  = [-5.0,0.4]; color_range['J-Ks']  = [0.07,0.80]
    FeH_range['BT-VT'] = [-2.7,0.4]; color_range['BT-VT'] = [0.19,1.49]
    FeH_range['VT-J']  = [-2.7,0.4]; color_range['VT-J']  = [0.77,2.56]
    FeH_range['VT-H']  = [-2.7,0.4]; color_range['VT-H']  = [0.77,3.16]
    FeH_range['VT-Ks'] = [-2.4,0.4]; color_range['VT-Ks'] = [0.99,3.29]
    FeH_range['b-y']   = [-3.7,0.5]; color_range['b-y']   = [0.18,0.72]

    coef = {}
    coef['B-V']   = [0.5665,0.4809,-0.0060,-0.0613,-0.0042,-0.0055]
    coef['V-RC']  = [0.4386,1.4614,-0.7014,-0.0807, 0.0142,-0.0015]
    coef['RC-IC'] = [0.3296,1.9716,-1.0225,-0.0298, 0.0329, 0.0035]
    coef['V-IC']  = [0.4033,0.8171,-0.1987,-0.0409, 0.0319, 0.0012]
    coef['V-J']   = [0.4669,0.3849,-0.0350,-0.0140, 0.0225, 0.0011]
    coef['V-H']   = [0.5251,0.2553,-0.0119,-0.0187, 0.0410, 0.0025]
    coef['V-Ks']  = [0.5057,0.2600,-0.0146,-0.0131, 0.0288, 0.0016]
    coef['J-Ks']  = [0.6393,0.6104, 0.0920,-0.0330, 0.0291, 0.0020]
    coef['BT-VT'] = [0.5839,0.4000,-0.0067,-0.0282,-0.0346,-0.0087]
    coef['VT-J']  = [0.4525,0.3797,-0.0357,-0.0082, 0.0123,-0.0009]
    coef['VT-H']  = [0.5286,0.2354,-0.0073,-0.0182, 0.0401, 0.0021]
    coef['VT-Ks'] = [0.4892,0.2634,-0.0165,-0.0121, 0.0249,-0.0001]
    coef['b-y']   = [0.5796,0.4812, 0.5747,-0.0633, 0.0042,-0.0055]

    extrapolation = kwargs.pop('extrapolation',False)

    if index not in coef:
        raise ColorIndexError(index, reference)

    if extrapolation or \
       (FeH_range[index][0]   <=  FeH  <= FeH_range[index][1] and \
        color_range[index][0] <= color <= color_range[index][1]):
        a = coef[index]
    else:
        raise ParamRangeError(index, color, reference)

    theta = a[0] + a[1]*color + a[2]*color**2 + a[3]*color*FeH \
            + a[4]*FeH + a[5]*FeH**2
    if index == 'b-y':
        M = [-1.9, 130.4, 125.7, 27.4]
        C = [-1003.7, 7325.9, -17207.4, 12977.7]
        for i in range(4):
            theta = theta +  M[i]*FeH**i + C[i]*color**i
    return 5040./theta
