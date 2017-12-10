import numpy.polynomial as poly

def temp_to_rgb(temp):
    '''
    Convert black-body temperature to RGB color.

    The colors are computed with piecewise polynomials. The coefficients are
    computed by fitting the CIE 1964 10-deg color matching functions given in
    `bbr_color.txt
    <http://www.vendian.org/mncharity/dir3/blackbody/UnstableURLs/bbr_color.html>`_
    on `Mitchell Charity's website
    <http://www.vendian.org/mncharity/dir3/blackbody/>`_.

    .. figure:: ../examples/true_colors/temp_rgb.png
       :alt: Blackbody temperature and real color
       :align: center
       :width: 600px
       :figwidth: 600px

       Relation of black-body temperatures v.s. RGB colors.

    Args:
        temp (float): Temperature in unit of K.
    Returns:
        tuple: A tuple cotaining (R, G, B) colors in range of [0.0, 1.0]
    See also:
        * :ref:`example_temp_color`
        

    '''
    # calculate red
    if temp <= 6500:
        red = 1.0
    else:
        t = min(temp, 40000)
        coeff = [2.68618299,  0.61788354, -0.40093348,
                 0.27789786, -0.25412922,  0.12710676]
        domain = (6600., 40000.)
        r2 = poly.Polynomial(coeff, domain=domain)
        red = min(1./r2(t), 1.0)

    # calculate green
    coeff_g1 = [0.5913838, 0.45727368, -0.10269557]
    domain_g1 = (1000., 6500.)
    g1 = poly.Polynomial(coeff_g1, domain=domain_g1)
    coeff_g2 = [1.85123116,  0.24567430, -0.14597665,
                0.09532498, -0.17510781,  0.12400952]
    domain_g2 = (6600., 40000.)
    g2 = poly.Polynomial(coeff_g2, domain=domain_g2)

    if temp <= 700:
        green = 0.0
    elif temp <= 6500:
        green = max(g1(temp), 0.0)
    elif temp >= 6600:
        t = min(temp, 40000)
        green = 1./g2(t)
    else:
        green = min(g1(temp), 1./g2(temp))

    # calculate blue
    if temp <= 1900:
        blue = 0.0
    elif temp < 6600:
        coeff_b1 = [0.42779584, 0.55646917, 0.06439607, -0.06693479]
        domain_b1 = (2000., 6500.)
        b1 = poly.Polynomial(coeff_b1, domain=domain_b1)
        blue = b1(temp)
        blue = max(blue, 0.0)
        blue = min(blue, 1.0)
    else:
        blue = 1.0

    return red, green, blue

def wavelength_to_rgb(wavelength):
    '''
    Convert wavelength to RGB color.

    See http://www.physics.sfasu.edu/astro/color/spectra.html

    Args:
        wavelength (float): Wavelength in Angstrom
    Returns:
        tuple: (R, G, B) color in 0~255.

    '''
    if 3800 <= wavelength <= 4400:
        r, g, b = (4400.-wavelength)/(4400.-3800.), 0, 1
    elif 4400 <= wavelength <= 4900:
        r, g, b = 0, (wavelength-4400.)/(4900.-4400.), 1
    elif 4900 <= wavelength <= 5100:
        r, g, b = 0, 1, (5100.-wavelength)/(5100.-4900.)
    elif 5100 <= wavelength <= 5800:
        r, g, b = (wavelength-5100.)/(5800.-5100.), 1, 0
    elif 5800 <= wavelength <= 6450:
        r, g, b = 1, (6450.-wavelength)/(6450.-5800.), 0
    elif 6450 <= wavelength <= 7800:
        r, g, b = 1, 0, 0
    elif wv > 7800:
        r, g, b = 1, 0, 0
    if wv < 3800:
        r, g, b = max((4400.-wavelength)/(4400.-3800.),0), 0, 1
    if wv > 7000:
        s = 0.3 + 0.7*(7800.-wavelength)/(7800.-7000.)
    elif wv < 4200:
        s = 0.3 + 0.7*(wavelength-3800.)/(4200.-3800.)
    else:
        s = 1.0

    if s < 0:
        s = 0.

    gamma = 0.8
    red   = int((s*r)**gamma*255)
    green = int((s*g)**gamma*255)
    blue  = int((s*b)**gamma*255)
    return (red, green, blue)
