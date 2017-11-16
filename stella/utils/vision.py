import numpy.polynomial as poly

def temp_to_rgb(temp):
    '''
    Convert black-body temperature to RGB color
    Parameters
    -----------
    temp


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
