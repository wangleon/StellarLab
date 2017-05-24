#!/usr/bin/env python

"""Astronomical and physical constants"""

# physical constants
c = 299792458.          # m/s,          light of speed, exactly
h = 6.62606957e-34      # kg m^2/s,     planck constant
e = 1.60217657e-19      # C,            elementary charge
G = 6.67384e-11         # m^3/(kg s),   gravitational constant
sigma = 5.670373e-8     # kg/(s^3 K^4), Stefan...Boltzmann constant
alpha = 7.2973525698e-3 #               fine-structure constant
m_e = 9.10938215e-31    # kg,           electron rest mass
m_n = 1.674927351e-27   # kg,           neutron mass
m_p = 1.672621777e-27   # kg,           proton mass

# distances
AU = 149597870700.      # m,            astronomical unit
pc = 3.08567758e16      # m,            parsec in meter
ly = 9460730472580800.  # m,            light-year, exactly
d_earth_sun = 14961875300. # m,         mean Earth-Sun distance

# celestial bodies      
R_sun   = 695508000.    # m,            solar radius
R_earth =   6378136.    # m,            earth equatorial radius
R_jup   =  71492000.    # m,            jupiter equatorial radius
M_sun   = 1.9891e30     # kg,           solar mass
M_earth = 5.97219e24    # kg,           earth mass
M_jup   = 1.89813e27    # kg,           jupiter mass
L_sun   = 3.846e26      # W,            solar luminosity (by GONG project)
L0      = 3.055e28      # W,            luminiosity at Mbol = 0.0 (IAU 1997)
Mbol_sun = 4.75         # mag,          bolometric magnitude of the Sun
Mv_sun   = 4.82         # mag,          absolute magnitude in V band of the Sun
BC_sun   = -0.07        # mag,          bolometric correction in V band of the Sun

# coordinations
# north galactic pole (NGP) in equatorial coordinate J2000.0
ALPHA_NGP = 192.859481  # RA (J2000.0) of NGP in degree, = 12[h]51[m]26.2755[s]
DELTA_NGP =  27.128251  # Dec (J2000.0) of NGP in degree, = 27[d]07[m]41.704[m]
L_NCP     = 122.931918  # l of NCP (North Celestial Pole) in degree

# north galactic pole (NGP) in equatorial coordinate B1950.0
# from Blauuw, Gum, Pawsey, Westerhout, 1960
ALPHA_NGP_B1950 = 192.25 # RA (B1950.0) of NGP in degree,  = 12[h]49[m]
DELTA_NGP_B1950 = 27.40  # Dec (B1950.0) of NGP in degree


# Geodetic coordinates IAU 1976                                                  
EARTH_EQ_RADIUS = 6.378140e6  # [m]                                              
EARTH_PO_RADIUS = 6.356755e6  # [m]                                              
EARTH_FLATTENING = 0.00335281
