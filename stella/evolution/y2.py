import os
import math
import numpy as np
import astropy.io.fits as fits

from .base import interpolate_data, interpolate_param

class _Y2(object):
    '''
    Yale-Yonsei stellar evolution track and isochrones.
    '''

    _data_path = '%s/evolution/Y2_tracks.fits'%os.getenv('STELLA_DATA')

    _mass_nodes = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5,
         1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,
         3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.5, 5.0]
    
    _z_nodes = [0.00001, 0.0001, 0.0004, 0.001, 0.004,
                0.007, 0.01, 0.02, 0.04, 0.06, 0.08]

    _alpha_nodes = [0.0, 0.3, 0.6]

    _missing_nodes = [(2.3, 0.0001, 0.0), (2.1, 0.001,  0.0),
                      (4.2, 0.02,   0.0), (4.2, 0.04,   0.0),
                      (4.2, 0.06,   0.0), (2.6, 0.08,   0.0),
                      (4.0, 0.08,   0.0), (4.2, 0.08,   0.0),
                      (5.0, 0.08,   0.0),
                      ]
    _wrong_nodes = [(2.4, 0.00001, 0.0), (2.9, 0.0001,  0.0),
                    (3.8, 0.0004,  0.0), (1.9, 0.0004,  0.0),
                    (2.3, 0.001,   0.0), (2.8, 0.001,   0.0),
                    (4.0, 0.001,   0.0), (2.7, 0.007,   0.0),
                    (2.9, 0.01,    0.0), (4.0, 0.02,    0.0),
                    ]

    _bad_nodes = {
            (0.08,    0.0): [2.6, 4.0, 4.2, 5.0],
            (0.06,    0.0): [4.2],
            (0.04,    0.0): [4.2],
            (0.02,    0.0): [4.0, 4.2],
            (0.01,    0.0): [2.9],
            (0.001,   0.0): [2.1, 2.3, 2.8, 4.0],
            (0.0004,  0.0): [1.9, 3.8],
            (0.0001,  0.0): [2.3],
            (0.00001, 0.0): [2.4],
            }

    _ngrid = 150

    def __init__(self):
        self._track_data     = {}
        self._isochrone_data = {}

    def get_track(self, mass, z, alpha=0.0, n=150):
        '''Get evolution track for given *M*, *Z*, and α ehancement.

        Args:
            mass (float): Stellar mass.
            z (float): Metal component.
            alpha (float): α ehancement.
            n (integer): Nunmber of points in evolution track.

        Returns:
            tuple: A tuple containing:
        '''

        if alpha in self._alpha_nodes:
            track = self._get_track_of_alpha(mass, z, alpha)
        else:
            track_lst = []
            for _alpha in self._alpha_nodes:
                track = self._get_track_of_alpha(mass, z, _alpha)
                track_lst.append(track)
            track = interpolate_param(track_lst, self._alpha_nodes, alpha)

        return track


    def _get_trackid(self, mass, z, alpha):
        '''Get Track ID.

        Args:
            mass (float): Stellar mass.
            z (float): Metal component.
            alpha (float): α ehancement.
        Returns:
            tuple: A tuple of three integers (`mass_id`, `z_id`, `alpha_id`) as
            the track ID.
        '''
        mass_id  = int(round(mass*10))
        z_id     = int(round(z*1e5))
        alpha_id = int(round(alpha*10))
        return (mass_id, z_id, alpha_id)

    def _load_all_tracks(self):
        '''Load all Y2 track data.'''
        data = fits.getdata(self._data_path)
        for row in data:
            trackid = self._get_trackid(row['mass'], row['z'], row['alpha'])
            logg = math.log10(row['mass']) + 4*(row['logTeff']
                        - math.log10(5777)) - row['logL'] + 4.44
            self._track_data[trackid] = (row['logTeff'], row['logL'],
                                         row['age'], logg)

    def _get_track_of_alpha(self, mass, z, alpha):
        '''Get evolution track for given *M*, *Z*, and alpha, of which alpha
        must be a value in grid nodes.

        '''
        if alpha not in self._alpha_nodes:
            print('Error: Alpha = %g not in alpha_nodes')
            raise ValueError

        if z in self._z_nodes:
            track = self._get_track_of_z(mass, z, alpha)
        else:
            track_lst = []
            iz = self._get_inodes(self._z_nodes, z)
            z_lst = self._z_nodes[iz:iz+4]
            for _z in z_lst:
                track = self._get_track_of_z(mass, _z, alpha)
                track_lst.append(track)
            track = interpolate_param(track_lst, np.log10(z_lst), math.log10(z))
        return track

    def _get_track_of_z(self, mass, z, alpha):
        '''Get evolution track for given *M*, *Z*, and alpha, of which *Z* must
        be a value in grid nodes.

        '''
        if z not in self._z_nodes:
            raise ValueError

        if (z, alpha) in self._bad_nodes:
            mass_nodes = [m for m in self._mass_nodes
                            if m not in self._bad_nodes[(z, alpha)]]
        elif abs(alpha - 0.6)<1e-3:
            # tracks of alpha=0.6, mass=4.2 are missed.
            mass_nodes = [m for m in self._mass_nodes
                            if abs(m-4.2)>1e-3]
        else:
            mass_nodes = self._mass_nodes

        if mass in mass_nodes:
            track = self._get_track_of_mass(mass, z, alpha)
        else:
            track_lst = []
            im = self._get_inodes(mass_nodes, mass)
            m_lst = mass_nodes[im:im+4]
            for _m in m_lst:
                track = self._get_track_of_mass(_m, z, alpha)
                track_lst.append(track)
            track = interpolate_param(track_lst, m_lst, mass)
        return track

    def _get_track_of_mass(self, mass, z, alpha):
        '''Get evolution track for given *M*, *Z*, and alpha, of which *M* must
        be a value in grid nodes.

        '''
        if mass not in self._mass_nodes:
            print('Error: M = %g not in mass_nodes')
            raise ValueError

        trackid = self._get_trackid(mass, z, alpha)
        if len(self._track_data) == 0:
            self._load_all_tracks()
        if trackid in self._track_data:
            track = self._track_data[trackid]
        else:
            print('missing Track: mass=%g, z=%g, alpha=%g'%(mass, z, alpha))
            raise ValueError
        return track

    def _get_inodes(self, nodes, value):
        '''Get nodes when using 4-points interpolation.

        '''
        l = len(nodes)
        i0 = np.searchsorted(nodes,value)
        if i0==0 or i0==1 or i0==2:
            return 0
        elif i0==l-2 or i0==l-1 or i0==l:
            return l-4
        else:
            return i0-2

Y2 = _Y2()
