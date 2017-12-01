import os
import math
import numpy as np
import astropy.io.fits as fits
from scipy.interpolate import splprep, splev

from ..parameter.metal import feh_to_z
from ..utils.interpolation import newton

_data_path  = '%s/evolution/Geneva_tracks.fits'%os.getenv('STELLA_DATA')
_z_nodes    = [0.001, 0.004, 0.008, 0.02, 0.04, 0.1]
_mass_nodes = [0.8, 0.9,  1.0, 1.25,  1.5,  1.7,  2.0,  2.5,  3.0,  4.0,  5.0,
               7.0, 9.0, 10.0, 12.0, 15.0, 20.0, 25.0, 40.0, 60.0, 85.0,120.0
              ]
_missing_nodes = [
    (0.001, 10.0),
    (0.004, 10.0),
    (0.008,  9.0),
    (0.02,  10.0),
    (0.04,  10.0),
    (0.1,   10.0), (0.1,   85.0), (0.1,  120.0),
    ]

def get_param_grid():
    '''
    Args:
    Returns:
        tuple: grid of parameter space (*Z*, *M*:sub:`0`) of Geneva tracks.
    '''
    param_grid = {}
    for z in _z_nodes:
        mass_lst = [m for m in _mass_nodes if (z,m) not in _missing_nodes]
        param_grid[z] = mass_lst
    return param_grid

def read_track(mass0, z):
    '''Read an evolution track in Geneva database.
    Args:
        mass0 (float): Initial mass
        z (float): Metal content
    Returns:
        tuple: lists of (log\ *T*:sub:`eff`, log\ *L*, age, *M*)
    '''
    f = fits.open(_data_path)
    data = f[1].data
    f.close()
    mask1 = (data['m0']==mass0)
    mask2 = (data['z']==z)
    mask = mask1*mask2
    row = data[mask][0]
    n = row['n']
    logTeff_lst = row['logTeff'][0:n]
    logL_lst    = row['logL'][0:n]
    age_lst     = row['age'][0:n]
    mass_lst    = row['mass'][0:n]
    print(z, mass0)
    return (logTeff_lst, logL_lst, age_lst, mass_lst)


def interpolate_track(track, n, k=1):
    '''Interpolate the evolution track
    Args:
        track (tuple): Input track (log\ *T*:sub:`eff`, log\ *L*, age, *M*)
        n (int): Number of interpolated points
        k (int, optinal): Degree of interpolated polynomial. Default is 1
    Returns:
        tuple: lists of (*M*, log\ *T*:sub:`eff`, log\ *L*, age)
    '''
    tck, u = splprep([track[0], track[1], track[2], track[3]], s=0, k=1)
    newx   = np.linspace(0, 1, n)
    newt    = splev(newx, tck)
    return (newt[0], newt[1], newt[2], newt[3])

def get_track(mass0, z, n=None):
    '''Get an evolution track for given (*M*:sub:`0`, *Z*) by interpolating
    the Geneva evolution track database.
    Args:
        mass0 (float): Initial mass
        z (float): Metal content
        n (int, optional): number of interpolated points
    Returns:
        tuple: lists of (*M*, log\ *T*:sub:`eff`, log\ *L*, age)
    '''
    ngrid = 51
    param_grid = get_param_grid()

    iz = _get_inodes(_z_nodes, z)

    inter3 = np.zeros((ngrid, 4))
    inter2 = np.zeros((ngrid, 4, 4))
    for jz in range(iz, iz+4):
        z_node = _z_nodes[jz]

        if mass0 in param_grid[z_node]:
            # if mass0 is a mass node
            track = read_track(mass0=mass0, z=z_node)
            if track[0].size != ngrid:
                track = interpolate_track(track, n=ngrid)
            # copy the parameters to inter2
            for k in range(4):
                inter2[:, k, jz-iz] = track[k]
                # k = 0, 1, 2, 3, which are the 4 parameters of track
                # jz - iz = 0, 1, 2, 3, which are the 4 data nodes
        else:
            # if mass0 is not a mass node
            # interpolate mass
            inter1 = np.zeros((ngrid, 4, 4))
            im = _get_inodes(param_grid[z_node], mass0)

            for jm in range(im, im+4):
                mass_node = param_grid[z_node][jm]

                track = read_track(mass0=mass_node, z=z_node)
                if track[0].size != ngrid:
                    track = interpolate_track(track, n=ngrid)
                for k in range(4):
                    inter1[:, k, jm-im] = track[k]
                    # k = 0, 1, 2, 3, which are the 4 parameters of track
                    # jm - im = 0, 1, 2, 3, which are the 4 data nodes

            for k1 in range(ngrid):
                for k2 in range(4):
                    # prepare for interpolation
                    vectx = param_grid[z_node][im:im+4]
                    vecty = inter1[k1, k2, :]
                    inter2[k1, k2, jz-iz] = newton(vectx, vecty, mass0)

    for k1 in range(ngrid):
        for k2 in range(4):
            # prepare for interpolation
            vectx = np.log10(_z_nodes[iz:iz+4])
            vecty = inter2[k1, k2, :]
            inter3[k1, k2] = newton(vectx, vecty, math.log10(z))

    logTeff_lst = inter3[:,0]
    logL_lst    = inter3[:,1]
    age_lst     = inter3[:,2]
    mass_lst    = inter3[:,3]

    track = (logTeff_lst, logL_lst, age_lst, mass_lst)

    if n is not None and n != ngrid:
        return interpolate_track(track, n=n)
    else:
        return track

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
