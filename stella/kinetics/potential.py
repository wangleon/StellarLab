import math
from ..constant import G, M_sun, pc

kpc = pc*1e3
varG = G*M_sun/kpc**2

class Potential(object):
    '''
    General class for potential.
    '''
    def __init__(self):
        pass

class PointPotential(Potential):
    '''
    Class for point potential
    '''

    def __init__(self, M):
        self.M = M

    def get_acce_cartesian(self, x, y, z):
        '''
        Get acceleration at given (*x*, *y*, *z*) in cartesian coordinates.

        Args:
            x (float): Galactic position *x* in unit of kpc
            y (float): Galactic position *y* in unit of kpc
            z (float): Galactic position *z* in unit of kpc
        Returns:
            tuple: acceleration components along (*x*, *y*, *z*) axes in unit of
                km s\ :sup:`−1`
        '''
        r = math.sqrt(x**2 + y**2 + z**2)
        coeff = varG*self.M/r**2
        ax = -coeff*x/r*1e-3
        ay = -coeff*y/r*1e-3
        az = -coeff*z/r*1e-3
        return (ax, ay, az)

    def get_vcirc(self, r):
        '''
        Get circular velocity at given distance *r*.

        Args:
            r (float): Distance to the Galactic center in unit of kpc
        Returns:
            float: circular velocity in unit of km s\ :sup:`−1`
        '''
        a_circ = varG*self.M/r**2
        return math.sqrt(a_circ*r*kpc)*1e-3

class HernquistPotential(Potential):
    '''
    Class for spherically symmetric Herquist potential.
    '''
    def __init__(self, M, a):
        self.M = M
        self.a = a

    def get_acce_cartesian(self, x, y, z):
        '''
        Get acceleration at given (*x*, *y*, *z*) in cartesian coordinates.

        Args:
            x (float): Galactic position *x* in unit of kpc
            y (float): Galactic position *y* in unit of kpc
            z (float): Galactic position *z* in unit of kpc
        Returns:
            tuple: acceleration components along (*x*, *y*, *z*) axes in unit of
                km s\ :sup:`−1`
        '''
        r = math.sqrt(x**2 + y**2 + z**2)
        coeff = varG*self.M/(r + self.a)**2
        ax = -coeff*x/r*1e-3
        ay = -coeff*y/r*1e-3
        az = -coeff*z/r*1e-3
        return (ax, ay, az)

    def get_vcirc(self, r):
        '''
        Get circular velocity at given distance *r*.

        Args:
            r (float): Distance to the Galactic center in unit of kpc
        Returns:
            float: circular velocity in unit of km s\ :sup:`−1`
        '''
        a_circ = varG*self.M/(r + self.a)**2
        return math.sqrt(a_circ*r*kpc)*1e-3

class MiyamotoNagaiPotential(Potential):
    '''
    Class for Miyamoto & Nagai potential.
    '''
    def __init__(self, M, a, b):
        self.M = M
        self.a = a
        self.b = b

    def get_acce_cartesian(self, x, y, z):
        '''
        Get acceleration at given (*x*, *y*, *z*) in cartesian coordinates.

        Args:
            x (float): Galactic position *x* in unit of kpc
            y (float): Galactic position *y* in unit of kpc
            z (float): Galactic position *z* in unit of kpc
        Returns:
            tuple: acceleration components along (*x*, *y*, *z*) axes in unit of
                km s\ :sup:`−1`
        '''
        zb = math.sqrt(z**2 + self.b**2)
        azb = (self.a + zb)**2
        xyazb = (x**2 + y**2 + azb)**1.5
        coeff = varG*self.M
        ax = -coeff*x/xyazb*1e-3
        ay = -coeff*y/xyazb*1e-3
        az = -coeff*z/zb*(self.a + zb)/xyazb*1e-3
        return (ax, ay, az)

    def get_vcirc(self, r):
        '''
        Get circular velocity at given distance *r* in the disk (*z* = 0).

        Args:
            r (float): Distance to the Galactic center in unit of kpc
        Returns:
            float: circular velocity in unit of km s\ :sup:`−1`
        '''
        a_circ = varG*self.M*r/(r**2 + (self.a + self.b)**2)**1.5
        return math.sqrt(a_circ*r*kpc)*1e-3

class NFWPotential(Potential):
    '''
    Class for spherically symmetric Navarro-Frenk-White potential.
    '''

    def __init__(self, M, rs):
        self.M  = M
        self.rs = rs

    def get_acce_cartesian(self, x, y, z):
        '''
        Get acceleration at given (*x*, *y*, *z*) in cartesian coordinates.

        Args:
            x (float): Galactic position *x* in unit of kpc
            y (float): Galactic position *y* in unit of kpc
            z (float): Galactic position *z* in unit of kpc
        Returns:
            tuple: acceleration components along (*x*, *y*, *z*) axes in unit of
                km s\ :sup:`−1`
        '''
        r = math.sqrt(x**2 + y**2 + z**2)
        coeff1 = varG*self.M
        coeff2 = 1./r/(r + self.rs) - 1./r**2*math.log(1. + r/self.rs)
        coeff = coeff1*coeff2
        ax = coeff*x/r*1e-3
        ay = coeff*y/r*1e-3
        az = coeff*z/r*1e-3
        return (ax, ay, az)

    def get_vcirc(self, r):
        '''
        Get circular velocity at given distance (*r*).

        Args:
            r (float): Distance to the Galactic center in unit of kpc
        Returns:
            float: circular velocity in unit of km s\ :sup:`−1`
        '''
        a_circ = -varG*self.M*(1./r/(r + self.rs) - 1./r**2*math.log(1. + r/self.rs))
        return math.sqrt(a_circ*r*kpc)*1e-3
