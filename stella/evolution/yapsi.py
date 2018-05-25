import os
import numpy as np
import astropy.io.fits as fits

from .base import interpolate_data, interpolate_param

class _YaPSI(object):
    '''YaPSI stellar evolution tracks.


    '''
    _y_nodes   = [0.25, 0.28, 0.31, 0.34, 0.37]
    _feh_nodes = [-1.5, -1.0, -0.5, 0.0, +0.3]
    _xz_nodes  = (
        ((0.749455,0.000545),(0.719477,0.000523),(0.689499,0.000501),
            (0.659520,0.000480),(0.629542,0.000458)), # [Fe/H] = -1.5
        ((0.748279,0.001721),(0.718348,0.001652),(0.688417,0.001583),
            (0.658485,0.001515),(0.628554,0.001446)), # [Fe/H] = -1.0
        ((0.744584,0.005416),(0.714801,0.005199),(0.685018,0.004982),
            (0.655234,0.004766),(0.625451,0.004549)), # [Fe/H] = -0.5
        ((0.733138,0.016862),(0.703812,0.016188),(0.674487,0.015513),
            (0.645161,0.014839),(0.615836,0.014164)), # [Fe/H] = +0.0
        ((0.717092,0.032908),(0.688408,0.031592),(0.659725,0.030275),
            (0.631041,0.028959),(0.602357,0.027643)), # [Fe/H] = +0.3
    )
    _mass_nodes = np.concatenate((
        np.arange(0.15, 0.40, 0.01),
        np.arange(0.40, 0.90, 0.02),
        np.arange(0.90, 1.80, 0.05),
        np.arange(1.80, 3.00, 0.10),
        np.arange(3.00, 5.01, 0.20),
    ))

    _alpha1, _alpha2 = 1.82126, 1.91804

    def __init__(self):
        self._track_data = {}

    def _load_tracks(self):
        '''Read evoution tracks.
        '''

    def _get_trackid(self, y, feh, mass, alpha):
        '''Get Track ID.

        Args:
            y (float): Initial helium content.
            feh (float): Metallcity.
            mass (float): Stellar mass.
            alpha (float): Mixing length.
        Returns:
            tuple: A tuple of four integers as the track ID.
        '''
        y_id    = int(round(y*100))
        feh_id  = int(round(feh*10))
        mass_id = int(round(mass*100))
        if abs(alpha-1.82126)<1e-5:
            alpha_id = 1
        elif abs(alpha-1.91804)<1e-5:
            alpha_id = 2
        else:
            print('Alpha ID does not exist')
            raise ValueError

        return (y_id, feh_id, mass_id, alpha_id)

    def get_track(self, y, feh, mass, n, minage=0):
        '''

        '''

        if y in self._y_nodes:
            track = self._get_track_of_feh(y, feh, mass, n, minage)
        else:
            track_lst = []
            iy = self._get_inodes(self._y_nodes, y)
            y_lst = self._y_nodes[iy:iy+4]
            for _y in y_lst:
                track = self._get_track_of_feh(_y, feh, mass, n, minage)
                track_lst.append(track)
            track = interpolate_param(track_lst, y_lst, y)

        return track
        
    def _get_inodes(self, nodes, value):
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

    def _has_alpha1(self, mass):
        '''
        '''
        return mass >= 0.60

    def _has_alpha2(self, ):
        '''
        '''
        return mass <= 1.10

    def _load_track(self, y, feh):
        '''
        '''
        data_path  = '%s/evolution/YaPSI'%os.getenv('STELLA_DATA')
        iy   = self._y_nodes.index(y)
        ifeh = self._feh_nodes.index(feh)
        x,z  = self._xz_nodes[ifeh][iy]
        filename = ('X%8.6f_Z%8.6f'%(x,z)).replace('.','p')+'.fits'
        filepath = os.path.join(data_path, filename)
        data = fits.getdata(filepath)
        for alpha in (self._alpha1, self._alpha2):
            mask1 = np.abs(data['alpha']-alpha)<1e-4
            for mass in self._mass_nodes:
                mask2 = np.abs(data['mass']-mass)<1e-4
                mask = mask1*mask2
                tracks = data[mask]
                trackid = self._get_trackid(y, feh, mass, alpha)
                track = (tracks['logTeff'], tracks['logL'], tracks['age'], tracks['logg'])
                self._track_data[trackid] = track

    def _get_track_of_mass(self, y, feh, mass, n, minage):
        if y not in self._y_nodes:
            print('Error: Y = %g not in y_nodes')
            raise ValueError
        if feh not in self._feh_nodes:
            print('Error: [Fe/H] = %g not in feh_nodes')
            raise ValueError

        if mass in self._mass_nodes:
            alpha = (self._alpha1, self._alpha2)[mass<=1.1]
            trackid = self._get_trackid(y, feh, mass, alpha)
            track = self._track_data[trackid]
            mask = track[2] > minage
            track = (track[0][mask], track[1][mask], track[2][mask], track[3][mask])
            track = interpolate_data(track, n)
        else:
            track_lst = []
            imass = self._get_inodes(self._mass_nodes, mass)
            mass_lst = self._mass_nodes[imass:imass+4]
            for _mass in mass_lst:
                alpha = (self._alpha1, self._alpha2)[_mass<=1.1]
                trackid = self._get_trackid(y, feh, _mass, alpha)
                track = self._track_data[trackid]
                mask = track[2] > minage
                track = (track[0][mask], track[1][mask], track[2][mask], track[3][mask])
                track = interpolate_data(track, n)
                track_lst.append(track)
            track = interpolate_param(track_lst, mass_lst, mass)
        return track

    def _get_track_of_feh(self, y, feh, mass, n, minage):
        if y not in self._y_nodes:
            print('Error: Y = %g not in y_nodes')
            raise ValueError

        if feh in self._feh_nodes:
            # load missing tracks
            trackid0 = self._get_trackid(y, feh, 1.0, self._alpha1)
            if trackid0 not in self._track_data:
                self._load_track(y, feh)
            track = self._get_track_of_mass(y, feh, mass, n, minage)
        else:
            track_lst = []
            ifeh = self._get_inodes(self._feh_nodes,feh)
            feh_lst = self._feh_nodes[ifeh:ifeh+4]
            for _feh in feh_lst:
                # load missing tracks
                trackid0 = self._get_trackid(y, _feh, 1.0, self._alpha1)
                if trackid0 not in self._track_data:
                    self._load_track(y, _feh)
                track = self._get_track_of_mass(y, _feh, mass, n, minage)
                track_lst.append(track)
            track = interpolate_param(track_lst, feh_lst, feh)
        return track

YaPSI = _YaPSI()
