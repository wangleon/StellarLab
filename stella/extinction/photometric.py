import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

def get_Stromgren_Eby(by,m1,c1,beta):
    """
    get E(b-y) based on the method of Olsen, 1988, A&A, 189, 173
    """

    # interpolated table from Crawford, 1975, AJ, 80, 955 table 1
    beta_lst = np.arange(2.590,2.720+1e-6, 0.01)
    m1_lst   = np.array([0.177,0.174,0.172,0.171,0.170,0.171,0.174,0.178,0.183,
                         0.189,0.196,0.204,0.214,0.226])
    c1_lst   = np.array([0.580,0.560,0.530,0.495,0.465,0.440,0.415,0.390,0.370,
                         0.350,0.330,0.310,0.290,0.270])
    f_m1 = intp.InterpolatedUnivariateSpline(beta_lst,m1_lst[::-1],k=3)
    f_c1 = intp.InterpolatedUnivariateSpline(beta_lst,c1_lst[::-1],k=3)
    m1_in = f_m1(beta)
    c1_in = f_c1(beta)

    m0,c0 = m1,c1

    delta_beta = 2.720 - beta

    niter = 0
    while(True):
        niter += 1
        delta_m0 = m1_in - m0
        delta_c0 = c0 - c1_in

        C = 4.9*delta_beta + 32.2*delta_m0 - 262.*delta_m0**2 - 1.31

        C = min(C, 1.6*delta_beta)
        if delta_m0 > 0.08:
            C = max(C, 0.13)
        else:
            C = max(C, -0.05)

        if delta_m0 < 0.06:
            D = (0.16 + 4.5*delta_m0 + 3.5*delta_beta)*delta_m0
        else:
            D = 0.24*delta_m0 + 0.035

        by0 = 0.217 + 1.34*delta_beta + 1.6*delta_beta**2 + C*delta_c0 - D

        Eby = by - by0
        m0_HJ = m1 + 0.30*Eby
        c0_HJ = c1 - 0.20*Eby

        if abs(m0_HJ-m0)<1e-5 and abs(c0_HJ-c0)<1e-5:
            break
        elif niter>=10:
            break
        else:
            m0,c0 = m0_HJ,c0_HJ

    return Eby


