from __future__ import print_function
import os
import math
import numpy as np
import astropy.io.fits as fits
from scipy.interpolate import splprep, splev

from ..parameter.metal import feh_to_z
from ..utils.interpolation import newton


_data_path = '%s/evolution/y2_tracks.fits'%os.getenv('STELLA_DATA')
_data = None

_mass_nodes = np.array([0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4,
    1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,
    3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.5, 5.0])

_z_nodes = np.array([0.00001, 0.0001, 0.0004, 0.001, 0.004, 0.007, 0.01, 0.02,
    0.04, 0.06, 0.08])

_alpha_nodes = np.array([0.0, 0.3, 0.6])

def _read_track(mass,z,alpha=0.0):
    '''
    Read Y2 evolution track for given alpha, mass and z
    '''
    if _data == None:
        _data = fits.getdata(_data_path)
    # search the track
    for row in _data:
        if alpha == row['alpha'] and mass == row['mass'] and z == row['z']:
            return row['logTeff'], row['logL'], row['age']
    print('Warning: the given track does not exist')
    return None, None, None

def _is_problematic_grid(mass,z,alpha):

    if alpha == 0.0:
        return (z,mass) in [(0.00001,2.4), (0.0001,2.9), (0.0004,3.8),
            (0.0004,1.9), (0.001,2.3), (0.001,2.8), (0.001,4.0), (0.007,2.7),
            (0.01,2.9), (0.02,4.0)]
    else:
        return False

def _in_grid(mass,z,alpha):
    if alpha not in _alpha_nodes or mass not in _mass_nodes or \
        z not in _z_nodes:
        return False

    if alpha == 0.0:
        return (z,mass) not in [(0.0001,2.3), (0.001,2.1), (0.02,4.2),
            (0.04,4.2), (0.06,4.2), (0.08,2.6), (0.08,4.0), (0.08,4.2),
            (0.08,5.0)]
    elif alpha == 0.6:
        return mass != 4.2
    else:
        return True


class Y2Track(object):

    _data_path = '%s/evolution/y2_tracks.fits'%os.getenv('STELLA_DATA')
    _mass_nodes = np.array([
        0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3,
        1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3,
        2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.4, 3.6,
        3.8, 4.0, 4.2, 4.5, 5.0])

    _z_nodes = np.array([
        0.00001, 0.0001, 0.0004, 0.001, 0.004,
        0.007, 0.01, 0.02, 0.04, 0.06, 0.08])

    _alpha_nodes = np.array([0.0,0.3,0.6])

    _ngrid = 150

    def __init__(self,**kwargs):
        '''args:
             mass  - stellar mass
             alpha - [alpha/Fe], default 0.0
             FeH   - [Fe/H]
             z     - metal component of star
             either FeH or z should be given
        '''

        #parse mass, [alpha/Fe], [Fe/H], z
        self.mass  = kwargs.pop('mass',None)
        self.alpha = kwargs.pop('alpha',0.0)
        self.feh   = kwargs.pop('feh', None)
        if self.feh == None:
            self.z = kwargs.pop('z', None)
            if self.z == None:
                raise ValueError
        else:
            self.z = feh_to_z(feh=self.feh, alpha=self.alpha)
        self.npoints = kwargs.pop('npoints',0)

        if None in [self.mass,self.z]:
            raise ValueError

        self._calc_track()

    def _read_track(self,mass,z,alpha=0.0):
        '''read evolution track for given alpha, mass and z'''

        hdulist = fits.open(self._data_path)
        data = hdulist[1].data
        hdulist.close()

        find = False
        for row in data:
            if abs(row['z']-z)<1e-7 and abs(row['mass']-mass)<1e-4:
                logTeff_lst = row['logTeff']
                logL_lst    = row['logL']
                age_lst     = row['age']
                find = True
                break

        if find:
            if age_lst.size != self._ngrid:
                tck, u = splprep([logTeff_lst,logL_lst,age_lst],s=0,k=3)
                newt   = np.linspace(0,1,self._ngrid)
                tmp    = splev(newt, tck)
                logTeff_lst = tmp[0]
                logL_lst    = tmp[1]
                age_lst     = tmp[2]
        else:
            logTeff_lst, logL_lst, age_lst = None, None, None

        return (logTeff_lst, logL_lst, age_lst)

    def _calc_track(self):

        # initialize mass list (m_lst)
        m_lst = [[] for s in self._z_nodes]
        # make m_lst
        for i,z in enumerate(self._z_nodes):
            mask = np.ones_like(self._mass_nodes)>0
            for j,mass in enumerate(self._mass_nodes):
                mask[j] = self._ingrid(mass,z,0.0) and \
                    not self._is_problematic_grid(mass,z,0.0)
            m_lst[i] = self._mass_nodes[mask]

        a_lst = self._alpha_nodes
        iz = self._get_inodes(self._z_nodes, self.z)

        inter3 = np.zeros((self._ngrid,3))
        inter2 = np.zeros((self._ngrid,3,4))
        for jz in np.arange(iz,iz+4):
            z_node = self._z_nodes[jz]

            if self.mass in m_lst[jz]:
                # if self.mass is a mass node
                track = self._read_track(mass=self.mass, z=z_node, alpha=0.0)
                for k in np.arange(3):
                    inter2[:,k,jz-iz] = track[k]
                    # k = 0, 1, 2. which is the 3 columns of track
                    # jz - iz = 0, 1, 2, 3. which is the 4 data nodes

            else:
                # if self.mass is not a mass node
                inter1 = np.zeros((self._ngrid,3,4))
                im = self._get_inodes(m_lst[jz],self.mass)

                for jm in range(im,im+4):
                    mass_node = m_lst[jz][jm]

                    track = self._read_track(mass=mass_node, z=z_node, alpha=0.0)
                    for k in np.arange(3):
                        inter1[:,k,jm-im] = track[k]
                        # k = 0, 1, 2.
                        # jm - im = 0, 1, 2, 3

                for k1 in range(self._ngrid):
                    for k2 in range(3):
                        # prepare for interpolation
                        vectx = m_lst[jz][im:im+4]
                        vecty = inter1[k1,k2,:]
                        x     = self.mass
                        inter2[k1,k2,jz-iz] = newton(vectx,vecty,x)

        for k1 in np.arange(self._ngrid):
            for k2 in np.arange(3):
                # prepare for interpolation
                vectx = np.log10(self._z_nodes[iz:iz+4])
                vecty = inter2[k1,k2,:]
                x     = math.log10(self.z)
                inter3[k1,k2] = newton(vectx,vecty,x)

        self.logTeff_lst = inter3[:,0]
        self.logL_lst    = inter3[:,1]
        self.age_lst     = inter3[:,2]

        if self.npoints >= 2:
            tck,u = splprep([
                self.age_lst,
                self.logTeff_lst,
                self.logL_lst], s=0, k=3)
            out = splev(np.linspace(0,1,self.npoints),tck)
            self.age_lst     = out[0]
            self.logTeff_lst = out[1]
            self.logL_lst    = out[2]

        self.Teff_lst = np.power(10,self.logTeff_lst)
        self.L_lst    = np.power(10,self.logL_lst)

        # calculate R, logg, and Mbol
        self.R_lst    = self.L_lst**0.5/(self.Teff_lst/5777.)**2
        self.logg_lst = math.log10(self.mass) - 2.*np.log10(self.R_lst) + 4.437
        self.Mbol_lst = 4.75 - 2.5*self.logL_lst


    def _get_inodes(self, nodes, value):
        '''get nodes when using 4-points interpolation'''
        l = len(nodes)
        i0 = np.searchsorted(nodes,value)
        if i0==0 or i0==1 or i0==2:
            return 0
        elif i0==l-2 or i0==l-1 or i0==l:
            return l-4
        else:
            return i0-2

    def _ingrid(self,mass,z,alpha):
        if alpha not in self._alpha_nodes or \
            mass not in self._mass_nodes  or \
              z not in self._z_nodes:
            return False

        para = (z, mass)
        if abs(alpha-0.0)<1e-3:
            missing_grids = [
                (0.0001,2.3),
                (0.001, 2.1),
                (0.02,  4.2),
                (0.04,  4.2),
                (0.06,  4.2),
                (0.08,  2.6), (0.08,  4.0), (0.08,  4.2), (0.08,  5.0),
                ]
            return (z,mass) not in missing_grids
        elif abs(alpha-0.6)<1e-3 and (mass-4.2)<1e-3:
            return False
        else:
            return True

    def _is_problematic_grid(self,mass,z,alpha):

        if abs(alpha-0.0)<1e-3:
            problematic_grids = [
                    (0.00001,2.4),
                    (0.0001,2.9),
                    (0.0004,3.8), (0.0004,1.9),
                    (0.001,2.3),(0.001,2.8),(0.001,4.0),
                    (0.007,2.7),
                    (0.01, 2.9),
                    (0.02, 4.0),
                ]
            return (z,mass) in problematic_grids
        else:
            return False

