import os
import math
import numpy as np
import astropy.io.fits as fits

from .base import interpolate_data, interpolate_param
from ..parameter.metal import feh_to_z

class _Geneva(object):
    '''Geneva stellar evolution tracks and isochrones.

    .. csv-table:: Structure of Geneva track file
        :header: Key, Type, Unit, Description
        :widths: 30, 30, 30, 100

        z,       float32,   ,            Metal content
        m0,      float32,   *M*:sub:`⊙`, Initial mass
        n,       integer16, ,            Number of data points
        logTeff, float32,   ,            log\ :sub:`10`\ (*T*:sub:`eff`)
        logL,    float32,   ,            log\ :sub:`10`\ (*L*\ /\ *L*:sub:`⊙`)
        age,     float32,   Gyr,         Age
        mass,    float32,   *M*:sub:`⊙`, Actual mass
    
    .. csv-table:: Structure of Geneva isocrhone file
        :header: Key, Type, Unit, Description
        :widths: 30, 30, 30, 100

        z,       float32,   ,            Metal content
        logage,  float32,   ,            log\ :sub:`10`\ (age)
        n,       integer16, ,            Number of data points
        m0,      float32,   *M*:sub:`⊙`, Initial mass
        mass,    float32,   *M*:sub:`⊙`, Actual mass
        logTeff, float32,   ,            log\ :sub:`10`\ (*T*:sub:`eff`)
        logL,    float32,   ,            log\ :sub:`10`\ (*L*\ /\ *L*:sub:`⊙`)

    '''

    def __init__(self):
        self._track_data = None
        self._isocrhone_data = None

    def _get_param_grid(self):
        '''Return a paramer grid that is available in the database.
        '''

    def _load_tracks(self):
        '''Read the whole Geneva track data file.
        '''
        data_path  = '%s/evolution/Geneva_tracks.fits'%os.getenv('STELLA_DATA')
        data = {}
        table = fits.getdata(data_path)
        for row in table:
            z  = row['z']
            m0 = row['m0']
            n  = row['n']
            logTeff_lst = row['logTeff'][0:n]
            logL_lst    = row['logL'][0:n]
            age_lst     = row['age'][0:n]
            mass_lst    = row['mass'][0:n]
            trackid = self._get_trackid(z, m0)
            data[trackid] = (logTeff_lst, logL_lst, age_lst, mass_lst)
        self._track_data = data

    def _load_isochrones(self):
        '''Read the whole Geneva isochrone file.
        '''
        data_path = '%s/evolution/Geneva_isochrones.fits'%os.getenv('STELLA_DATA')
        data = {}
        table = fits.getdata(data_path)
        for row in table:
            z      = row['z']
            logage = row['logage']
            n      = row['n']
            m0_lst      = row['m0'][0:n]
            mass_lst    = row['mass'][0:n]
            logTeff_lst = row['logTeff'][0:n]
            logL_lst    = row['logL'][0:n]
            isochroneid = self._get_isochroneid(z, m0)
            data[isochroneid] = (m0_lst, mass_lst, logTeff_lst, logL_lst)
        self._isochrone_data = data

    def _get_trackid(self, z, mass0):
        '''Get Track ID.

        Args:
            mass0 (float): Initial mass
            z (float): Metal content
        Returns:
            tuple: Track ID
        '''
        return (int(round(z*1000)), int(round(mass0*100)))

    def _get_isochroneid(self, z, logage):
        '''Get Isochrone ID.

        Args:
            z (float): Metal content
            logage (float): log\ :sub:`10`\ (age)
        Returns:
            tuple: Isochrone ID.
        '''
        return (int(round(z*1000)), int(round(logage*100)))

    def get_track(self, mass0, z, n=None):
        '''Get an evolution track for given (*M*:sub:`0`, *Z*) by interpolating
        the Geneva evolution track database.
    
        Args:
            mass0 (float): Initial mass.
            z (float): Metal content.
            n (int, optional): number of interpolated points.
        Returns:
            tuple: A tuple containing (log\ *T*:sub:`eff`, log\ *L*, age, *M*).
        '''
        if self._track_data is None:
            self._load_tracks()

        z_nodes = [0.001, 0.004, 0.008, 0.02, 0.04, 0.1]
        m_nodes = [0.8, 0.9, 1.0, 1.25, 1.5, 1.7, 2.0, 2.5, 3.0, 4.0, 5.0,
                   7.0, 9.0, 10.0, 12.0, 15.0, 20.0, 25.0, 40.0, 60.0, 85.0,
                   120.0]
        missing_nodes = [(0.001, 10.0), (0.004, 10.0), (0.008, 9.0),
                         (0.02, 10.0), (0.04, 10.0), (0.1, 10.0), (0.1, 85.0),
                         (0.1, 120.0)]

        # get param grid
        param_grid = {}
        for _z in z_nodes:
            param_grid[_z] = [_m for _m in m_nodes
                              if (_z, _m) not in missing_nodes]

        ngrid = 51
    
        if z in param_grid:
            # input z in parameter grid
            if mass0 in param_grid[z]:
                # input mass0 in parameter grid
                trackid = self._get_trackid(z, mass0)
                track = self._track_data[trackid]
                if track[0].size != ngrid:
                    track = interpolate_data(track, n=ngrid)
            else:
                # input mass0 NOT in parameter grid. Interpolate over mass0 space
                im = _get_inodes(param_grid[z], mass0)
                mass0_lst = param_grid[z][im:im+4]
                track_lst = []
                for _mass0 in mass0_lst:
                    trackid = self._get_trackid(z, _mass0)
                    track = self._track_data[trackid]
                    if track[0].size != ngrid:
                        track = interpolate_data(track, n=ngrid)
                    track_lst.append(track)
                track = interpolate_param(track_lst, mass0_lst, mass0)
        else:
            # input z Not in parameter grid. Interpolate over log10(z) space
            iz = _get_inodes(z_nodes, z)
            z_lst = z_nodes[iz:iz+4]
            trackz_lst = []
            for _z in z_lst:
                if mass0 in param_grid[_z]:
                    # input mass0 in parameter grid
                    trackid = self._get_trackid(_z, mass0)
                    track = self._track_data[trackid]
                    if track[0].size != ngrid:
                        track = interpolate_data(track, n=ngrid)
                else:
                    # input mass0 NOT in parameter grid. Interpolate over mass0
                    # space
                    im = _get_inodes(param_grid[_z], mass0)
                    mass0_lst = param_grid[_z][im:im+4]
                    trackm_lst = []
                    for _mass0 in mass0_lst:
                        trackid = self._get_trackid(_z, _mass0)
                        track = self._track_data[trackid]
                        if track[0].size != ngrid:
                            track = interpolate_data(track, n=ngrid)
                        trackm_lst.append(track)
                    track = interpolate_param(trackm_lst, mass0_lst, mass0)
                trackz_lst.append(track)
            track = interpolate_param(trackz_lst, np.log10(z_lst), math.log10(z))
    
        if n is not None and n != ngrid:
            # interpolate for given number of points
            return interpolate_data(track, n=n)
        else:
            return track

    def get_isochrone(self, z, logage, n=None):
        '''Get an isochrone for given (*Z*, age) by interpolating the Geneva
        evolution isochrone database.

        Args:
            z (float): Metal content
            logage (float): log\ :sub:`10`\ (age)
            n (int, optional): number of interpolated points
        Returns:
            tuple: A tuple containing (*M*:sub:`0`, *M*, log\ *T*:sub:`eff`, log\ *L*)
        '''

        ngrid = 600

        if self._isochrone_data is None:
            self.load_isochrones()

        z_nodes = [0.001, 0.004, 0.008, 0.02, 0.04, 0.1]
        a_nodes = [3.00, 5.00, 5.30, 5.59, 5.80, 5.90, 6.00, 6.05, 6.09, 6.15,
                   6.19, 6.25, 6.30, 6.34, 6.40, 6.44, 6.50, 6.55, 6.59, 6.65,
                   6.69, 6.75, 6.80, 6.84, 6.90, 6.94, 7.00, 7.05, 7.09, 7.15,
                   7.19, 7.25, 7.30, 7.34, 7.40, 7.44, 7.50, 7.55, 7.59, 7.65,
                   7.69, 7.75, 7.80, 7.84, 7.90, 7.94, 8.00, 8.05, 8.10, 8.14,
                   8.19, 8.25, 8.30, 8.35, 8.39, 8.44, 8.50, 8.55, 8.60, 8.64,
                   8.69, 8.75, 8.80, 8.85, 8.89, 8.94, 9.00, 9.05, 9.10, 9.14,
                   9.19, 9.25, 9.30, 9.35, 9.39, 9.44, 9.50, 9.55, 9.60, 9.64,
                   9.69, 9.75, 9.80, 9.85, 9.89, 9.94, 10.00, 10.05, 10.10,
                   10.14, 10.19]

        if z in z_nodes:
            if logage in a_nodes:
                isochroneid = self._get_isochroneid(z, logage)
                isochrone = self._isochrone_data[isochroneid]
                if isochrone[0].size != ngrid:
                    isochrone = interpolate_data(isochrone, n=ngrid)
            else:
                ia = _get_inodes(a_nodes, logage)
                logage_lst = a_nodes[ia:ia+4]
                isochrone_lst = []
                for _logage in logage_lst:
                    isochroneid = self._get_isochroneid(z, _logage)
                    isochrone = self._isochrone_data[isochroneid]
                    if isochrone[0].size != ngrid:
                        isochrone = interpolate_data(isochrone, n=ngrid)
                    isochrone_lst.append(isochrone)
                isochrone = interpolate_param(isochrone_lst, logage_lst, logage)
        else:
            iz = _get_inodes(z_nodes, z)
            z_lst = z_nodes[iz:iz+4]
            isochronez_lst = []
            for _z in z_lst:
                if logage in a_nodes:
                    isochroneid = self._get_isochroneid(_z, logage)
                    isochrone = self._isochrone_data[isochroneid]
                    if isochrone[0].size != ngrid:
                        isochrone = interpolate_data(isochrone, n=ngrid)
                else:
                    ia = _get_inodes(a_nodes, logage)
                    logage_lst = a_nodes[ia:ia+4]
                    isochroneage_lst = []
                    for _logage in logage_lst:
                        isochroneid = self._get_isochroneid(_z, _logage)
                        isochrone = self._isochrone_data[isochroneid]
                        if isochrone[0].size != ngrid:
                            isochrone = interpolate_data(isochrone, n=ngrid)
                        isochroneage_lst.append(isochrone)
                    isochronez = interpolate_param(isochroneage_lst, logage_lst, logage)
                isochronez_lst.append(isochronez)
            isochrone = interepolate_param(isochronez_lst, np.log10(z_lst), math.log10(z))

        return isochrone

Geneva = _Geneva()


def _get_inodes(nodes, value):
    '''Get the begining index of the 4-points interpolation.

    Args:
        nodes (list): Input node list
        value (intger or float): Input value
    Returns:
        int: Beginning index of 4-points interpolation
    '''
    i0 = np.searchsorted(nodes, value)
    i = i0-2
    i = max(i, 0)
    i = min(i, len(nodes)-4)
    return i

