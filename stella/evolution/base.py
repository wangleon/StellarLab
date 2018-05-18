import numpy as np
from scipy.interpolate import splprep, splev
from ..utils.interpolation import newton

def interpolate_data(track, n, k=1):
    '''Interpolate the evolution track.

    Args:
        track (tuple): Input track as a tuple of numpy arrays.
        n (integer): Number of interpolated points
        k (integer): Degree of interpolated polynomial. Default is 1
    Returns:
        tuple: A tuple of numpy arrays.
    '''
    tck, u = splprep([track[0], track[1], track[2], track[3]], s=0, k=1)
    newx   = np.linspace(0, 1, n)
    newt   = splev(newx, tck)
    return (newt[0], newt[1], newt[2], newt[3])


def interpolate_param(track_lst, param_lst, param):
    '''Interpolate the tracks over a certain parameter space.

    Args:
        track_lst (list): List of track tuples
        param_lst (list): List of node parameters in grid
        param (integer or float): Input parameter
    Returns:
        tuple: A tuple containing (log\ *T*:sub:`eff`, log\ *L*, age, *M*)
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

