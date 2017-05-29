import math
from ..constant import G, M_sun, pc

kpc = pc*1e3
rG = G*M_sun/kpc**2

class Potential(object):

    def __init__(self):
        pass

class PointPotential(Potential):

    def __init__(self, M):
        self.M = M

    def acce_cartesian(self, x, y, z):
        r = math.sqrt(x**2 + y**2 + z**2)
        ax = -rG*self.M/r**3*x
        ay = -rG*self.M/r**3*y
        az = -rG*self.M/r**3*z
        return (ax, ay, az)

class HernquistPotential(Potential):

    def __init__(self, M, a):
        self.M = M
        self.a = a

    def acce_cartesian(self, x, y, z):
        r = math.sqrt(x**2 + y**2 + z**2)
        ax = -rG*self.M/(r + self.a)**2*x/r
        ay = -rG*self.M/(r + self.a)**2*y/r
        az = -rG*self.M/(r + self.a)**2*z/r
        return (ax, ay, az)

class MiyamotoNagaiPotential(Potential):

    def __init__(self, M, a, b):
        self.M = M
        self.a = a
        self.b = b

    def acce_cartesian(self, x, y, z):
        zb = math.sqrt(z**2 + self.b**2)
        azb = (self.a + zb)**2
        xyazb = (x**2 + y**2 + azb)**1.5
        ax = -rG*self.M*x/xyazb
        ay = -rG*self.M*y/xyazb
        az = -rG*self.M*z/zb*(self.a + zb)/xyazb
        return (ax, ay, az)

class NFWPotential(Potential):

    def __init__(self, M, rs):
        self.M  = M
        self.rs = rs

    def acce_cartesian(self, x, y, z):
        r = math.sqrt(x**2 + y**2 + z**2)
        coeff = 1./r/(r + self.rs) - 1./r**2*math.log(1. + r/self.rs)
        ax = rG*self.M*coeff*x/r
        ay = rG*self.M*coeff*y/r
        az = rG*self.M*coeff*z/r
        return (ax, ay, az)
