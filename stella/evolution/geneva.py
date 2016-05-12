#!/usr/bin/env python
import os
import math
import numpy as np
import astropy.io.fits as fits
from scipy.interpolate import splprep, splev

from ..param.metal import feh_to_z
from ..utils.interpolation import newton

class GenevaTrack(object):
    _data_path = '%s/evolution/geneva_tracks.fits'%os.getenv('STELLA_DATA')
    _mass_nodes = np.array([
        0.8, 0.9,  1.0, 1.25,  1.5,  1.7,  2.0,  2.5,  3.0,  4.0,  5.0,
        7.0, 9.0, 10.0, 12.0, 15.0, 20.0, 25.0, 40.0, 60.0, 85.0,120.0])

    _z_nodes  = np.array([0.001,0.004,0.008,0.02,0.04,0.1])

    _ngrid = 51

    def __init__(self,**kwargs):
        '''args:
                mass   - stellar mass
                feh    - [Fe/H]
                z      - metal component of star
                either FeH or z should be given
        '''

        #parse mass, [alpha/Fe], [Fe/H], z
        self.mass0 = kwargs.pop('mass',None)
        self.alpha = kwargs.pop('alpha',0.0)
        self.feh   = kwargs.pop('feh', None)
        if self.feh == None:
            self.z = kwargs.pop('z', None)
            if self.z == None:
                raise ValueError
        else:
            self.z = feh_to_z(feh=self.feh,  alpha=self.alpha)
        self.npoints = kwargs.pop('npoints',0)

        self._calc_track()

    def _calc_track(self):
        # initialize mass list (m_lst)
        m_lst = [[] for s in self._z_nodes]
        # make m_lst
        for i,z in enumerate(self._z_nodes):
            mask = np.ones_like(self._mass_nodes)>0
            for j,mass in enumerate(self._mass_nodes):
                mask[j] = self._ingrid(mass=mass, z=z)
            m_lst[i] = self._mass_nodes[mask]

        iz = self._get_inodes(self._z_nodes, self.z)

        inter3 = np.zeros((self._ngrid,4))
        inter2 = np.zeros((self._ngrid,4,4))
        for jz in range(iz,iz+4):
            z_node = self._z_nodes[jz]

            if self.mass0 in m_lst[jz]:
                # if self.mass0 is a mass node
                track = self._read_track(mass=self.mass0, z=z_node)
                for k in range(4):
                    inter2[:,k,jz-iz] = track[k]
                    # k = 0, 1, 2. which is the 3 columns of track
                    # jz - iz = 0, 1, 2, 3. which is the 4 data nodes

            else:
                # if self.mass0 is not a mass node
                inter1 = np.zeros((self._ngrid,4,4))
                im = self._get_inodes(m_lst[jz], self.mass0)

                for jm in range(im,im+4):
                    mass_node = m_lst[jz][jm]

                    track = self._read_track(mass=mass_node, z=z_node)
                    for k in range(4):
                        inter1[:,k,jm-im] = track[k]
                        # k = 0, 1, 2.
                        # jm - im = 0, 1, 2, 3

                for k1 in range(self._ngrid):
                    for k2 in range(4):
                        # prepare for interpolation
                        vectx = m_lst[jz][im:im+4]
                        vecty = inter1[k1,k2,:]
                        x     = self.mass0
                        inter2[k1,k2,jz-iz] = newton(vectx,vecty,x)

        for k1 in range(self._ngrid):
            for k2 in range(4):
                # prepare for interpolation
                vectx = np.log10(self._z_nodes[iz:iz+4])
                vecty = inter2[k1,k2,:]
                x     = math.log10(self.z)
                inter3[k1,k2] = newton(vectx,vecty,x)

        self.logTeff_lst = inter3[:,0]
        self.logL_lst    = inter3[:,1]
        self.age_lst     = inter3[:,2]
        self.mass_lst    = inter3[:,3]

        if self.npoints >= 2:
            tck,u = splprep([
                self.logTeff_lst, self.logL_lst, self.age_lst, self.mass_lst],
                u=np.arange(self._ngrid), s=0, k=1)
            out = splev(np.linspace(0,self._ngrid-1,self.npoints),tck)
            self.logTeff_lst = out[0]
            self.logL_lst    = out[1]
            self.age_lst     = out[2]
            self.mass_lst    = out[3]

        self.Teff_lst = np.power(10,self.logTeff_lst)
        self.L_lst    = np.power(10,self.logL_lst)
        
        # calculate R, logg, and Mbol
        self.R_lst    = self.L_lst**0.5/(self.Teff_lst/5777.)**2
        self.logg_lst = np.log10(self.mass_lst) - 2.*np.log10(self.R_lst) + 4.437
        self.Mbol_lst = 4.75 - 2.5*self.logL_lst

    def _ingrid(self,mass,z):
        missing_grids = [
            (0.001, 10.0), (0.004, 10.0), (0.008, 9.0),
            (0.02,  10.0), (0.04,  10.0),
            (0.1,  120.0), (0.1,   85.0), (0.1,  10.0)
            ]
        return (z,mass) not in missing_grids

    def _read_track(self,mass,z):

        hdulist = fits.open(self._data_path)
        data = hdulist[1].data

        find = False
        for row in data:
            if abs(row['z']-z)<1e-7 and abs(row['m0']-mass)<1e-4:
                n = row['n']
                logTeff_lst = row['logTeff'][0:n]
                logL_lst    = row['logL'][0:n]
                age_lst     = row['age'][0:n]
                mass_lst    = row['mass'][0:n]
                find = True
                break

        if find:
            if age_lst.size != self._ngrid:
                tck, u = splprep(
                        [logTeff_lst,logL_lst,age_lst,mass_lst],s=0,k=3)
                newt   = np.linspace(0,1,self._ngrid)
                tmp    = splev(newt, tck)
                logTeff_lst = tmp[0]
                logL_lst    = tmp[1]
                age_lst     = tmp[2]
                mass_lst    = tmp[3]
        else:
            logTeff_lst, logL_lst, age_lst = None, None, None, None

        return (logTeff_lst, logL_lst, age_lst, mass_lst)

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

