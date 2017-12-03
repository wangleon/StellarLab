import os
import math
import numpy as np
import astropy.io.fits as fits
from scipy.interpolate import splprep, splev

from ..parameter.metal import feh_to_z
from ..utils.interpolation import newton

class Geneva(object):
    '''Geneva stellar evolution tracks and isochrones.


    .. csv-table:: Structure of Geneva track file
        :header: "Key", "Type", "Unit", "Description"
        :widths: 30, 30, 30, 100

        "z",       "float32",   "",            "Metal content"
        "m0",      "float32",   "*M*:sub:`⊙`", "Initial mass"
        "n",       "integer16", "",            "Number of data points"
        "logTeff", "float32",   "",            "log\ :sub:`10`\ (*T*:sub:`eff`)"
        "logL",    "float32",   "",            "log\ :sub:`10`\ (*L*\ /\ *L*:sub:`⊙`)"
        "age",     "float32",   "Gyr",         "Age"
        "mass",    "float32",   "*M*:sub:`⊙`", "Actual mass"
    
    .. csv-table:: Structure of Geneva isocrhone file
        :header: "Key", "Type", "Unit", "Description"
        :widths: 30, 30, 30, 100

        "z",       "float32",   "",            "Metal content"
        "logage",  "float32",   "",            "log\ :sub:`10`\ (age)"
        "n",       "integer16", "",            "Number of data points"
        "m0",      "float32",   "*M*:sub:`⊙`", "Initial mass"
        "mass",    "float32",   "*M*:sub:`⊙`", "Actual mass"
        "logTeff", "float32",   "",            "log\ :sub:`10`\ (*T*:sub:`eff`)"
        "logL",    "float32",   "",            "log\ :sub:`10`\ (*L*\ /\ *L*:sub:`⊙`)"

    '''
    _trackdata_path  = '%s/evolution/Geneva_tracks.fits'%os.getenv('STELLA_DATA')

    def __init__(self):
        self._track_data = None

    def _get_param_grid(self):
        '''Return a paramer grid that is available in the database.
        '''

    def _load_tracks(self):
        '''Read the whole Geneva track data file.
        '''
        data = {}
        table = fits.getdata(self._trackdata_path)
        for row in table:
            z = row['z']
            m0 = row['m0']
            n = row['n']
            logTeff_lst = row['logTeff'][0:n]
            logL_lst    = row['logL'][0:n]
            age_lst     = row['age'][0:n]
            mass_lst    = row['mass'][0:n]
            trackid = self._get_trackid(z, m0)
            data[trackid] = (logTeff_lst, logL_lst, age_lst, mass_lst)
        self._track_data = data

    def _get_trackid(self, z, mass0):
        '''Get Track ID.

        Args:
            mass0 (float): Initial mass
            z (float): Metal content
        Returns:
            tuple: Track ID
        '''
        return (int(round(z*1000)), int(round(mass0*100)))

    def get_track(self, mass0, z, n=None):
        '''Get an evolution track for given (*M*:sub:`0`, *Z*) by interpolating
        the Geneva evolution track database.
    
        Args:
            mass0 (float): Initial mass
            z (float): Metal content
            n (int, optional): number of interpolated points
        Returns:
            tuple: lists of (log\ *T*:sub:`eff`, log\ *L*, age, *M*)
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
                    track = interpolate_track(track, n=ngrid)
            else:
                # input mass0 NOT in parameter grid. Interpolate over mass0 space
                im = _get_inodes(param_grid[z], mass0)
                mass0_lst = param_grid[z][im:im+4]
                track_lst = []
                for _mass0 in mass0_lst:
                    trackid = self._get_trackid(z, _mass0)
                    track = self._track_data[trackid]
                    if track[0].size != ngrid:
                        track = interpolate_track(track, n=ngrid)
                    track_lst.append(track)
                track = interpolate_param(track_lst, mass0_lst, mass0)
        else:
            # input z Not in parameter grid. Interpolate over log10(z) space
            iz = _get_inodes(self.z_nodes, z)
            z_lst = self.z_nodes[iz:iz+4]
            trackz_lst = []
            for _z in z_lst:
                if mass0 in param_grid[_z]:
                    # input mass0 in parameter grid
                    trackid = self._get_trackid(_z, mass0)
                    track = self._track_data[trackid]
                    if track[0].size != ngrid:
                        track = interpolate_track(track, n=ngrid)
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
                            track = interpolate_track(track, n=ngrid)
                        trackm_lst.append(track)
                    track = interpolate_param(trackm_lst, mass0_lst, mass0)
                trackz_lst.append(track)
            track = interpolate_param(trackz_lst, np.log10(z_lst), math.log10(z))
    
        if n is not None and n != ngrid:
            # interpolate for given number of points
            return interpolate_track(track, n=n)
        else:
            return track


def interpolate_track(track, n, k=1):
    '''Interpolate the evolution track.

    Args:
        track (tuple): Input track (log\ *T*:sub:`eff`, log\ *L*, age, *M*)
        n (int): Number of interpolated points
        k (int, optinal): Degree of interpolated polynomial. Default is 1
    Returns:
        tuple: lists of (log\ *T*:sub:`eff`, log\ *L*, age, *M*)
    '''
    tck, u = splprep([track[0], track[1], track[2], track[3]], s=0, k=1)
    newx   = np.linspace(0, 1, n)
    newt   = splev(newx, tck)
    return (newt[0], newt[1], newt[2], newt[3])


def _get_inodes(nodes, value):
    '''Get the begining index of the 4-points interpolation.

    Args:
        nodes (list): Input node list
        value (int or float): Input value
    Returns:
        int: Beginning index of 4-points interpolation
    '''
    i0 = np.searchsorted(nodes, value)
    i = i0-2
    i = max(i, 0)
    i = min(i, len(nodes)-4)
    return i

def interpolate_param(track_lst, param_lst, param):
    '''Interpolate the tracks over a certain parameter space.

    Args:
        track_lst (list): List of track tuples.
        param_lst (list): List of node parameters in grid.
        param (int or float): Input parameter
    Returns:
        tuple: lists of (log\ *T*:sub:`eff`, log\ *L*, age, *M*)
    '''

    ntrack = len(track_lst)
    nparam = len(track_lst[0])
    ngrid  = track_lst[0][0].size

    inter1 = np.zeros((ngrid, nparam, ntrack))
    inter2 = np.zeros((ngrid, nparam))

    for it, track in enumerate(track_lst):
        for ip, v_lst in enumerate(track):
            inter1[:, ip, it] = v_lst

    for k1 in range(ngrid):
        for k2 in range(nparam):
            inter2[k1, k2] = newton(param_lst, inter1[k1, k2, :], param)

    newtrack = tuple(inter2[:,k] for k in range(nparam))

    return newtrack

