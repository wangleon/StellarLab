import os
import numpy as np
import astropy.io.fits as fits

from .base import interpolate_data, interpolate_param
from ..utils.download import get_file

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

    # alm: alpha mixing length
    _alm1, _alm2 = 1.82126, 1.91804

    def __init__(self):
        self._track_data = {}

    def _load_tracks(self):
        '''Read evoution tracks.
        '''

    def _get_trackid(self, y, feh, mass, alm):
        '''Get Track ID.

        Args:
            y (float): Initial helium content.
            feh (float): Metallcity.
            mass (float): Stellar mass.
            alm (float): Mixing length.
        Returns:
            tuple: A tuple of four integers (`y_id`, `feh_id`, `mass_id`,
                `alpha_id`) as the track ID.
        '''
        y_id    = int(round(y*100))
        feh_id  = int(round(feh*10))
        mass_id = int(round(mass*100))
        if abs(alm-self._alm1)<1e-5:
            alm_id = 1
        elif abs(alm-self._alm2)<1e-5:
            alm_id = 2
        else:
            print('Alpha ID does not exist')
            raise ValueError

        return (y_id, feh_id, mass_id, alm_id)

    def get_track(self, y, feh, mass, n=0, minage=0):
        '''
        Get the evolution track for given *Y*, [Fe/H] and *M*.

        Args:
            y (float):
            feh (float): Metallicity.
            mass (float): Stellar mass.
            n (integer): Number of points in evolution track.
            minage (float): Minimum age returned in the evolution track.

        Returns:
            tuple: A tuple containing:

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
            nodes (list): Input node list.
            value (intger or float): Input value.
        Returns:
            integer: Beginning index of 4-points interpolation.
        '''
        i0 = np.searchsorted(nodes, value)
        i = i0-2
        i = max(i, 0)
        i = min(i, len(nodes)-4)
        return i

    def _load_track(self, y, feh, mass, alm):
        '''
        Load the track of given *y* and *feh* values and pack it into the cache.

        Args:
            y (float): Helium content (*Y*).
            feh (float): Metallicity ([Fe/H]).
            mass (float): Stellar mass (*M*).
            alm (float): Mixing-length (*alpha*).

        Notes:
            The track data is a tuple containing four arrays
            (log\ *T*:sub:`eff`, log\ *L*, age, log\ *g*).

        '''
        iy   = self._y_nodes.index(y)
        ifeh = self._feh_nodes.index(feh)
        x,z  = self._xz_nodes[ifeh][iy]

        data_path  = 'thirdpartydata/yapsi'
        folder = 'X{:8.6f}_Z{:8.6f}'.format(x, z).replace('.','p')
        fname = 'M{:4.2f}_X{:8.6f}_Z{:8.6f}_A{:7.5f}'.format(
                mass, x, z, alm).replace('.', 'p')+'.trk'
        filepath = os.path.join(data_path, folder, fname)
        filename = get_file(filepath)

        # load data
        data = []
        file1 = open(filename)
        for row in file1:
            if row[0]=='#':
                continue
            col = row.split()
            age     = float(col[2])
            logL    = float(col[6])
            logR    = float(col[7])
            logg    = float(col[8])
            logTeff = float(col[9])
            item = (logTeff, logL, age, logR, logg)
            data.append(item)
        file1.close()

        # repack parameters
        logTeff_lst = np.array([row[0] for row in data])
        logL_lst    = np.array([row[1] for row in data])
        age_lst     = np.array([row[2] for row in data])
        logR_lst    = np.array([row[3] for row in data])
        logg_lst    = np.array([row[4] for row in data])

        track = (logTeff_lst, logL_lst, age_lst, logg_lst, logR_lst)
        return track

    def _get_track_of_mass(self, y, feh, mass, n=0, minage=0):
        if y not in self._y_nodes:
            print('Error: Y = {} not in y_nodes'.format(y))
            raise ValueError
        if feh not in self._feh_nodes:
            print('Error: [Fe/H] = {} not in feh_nodes'.format(feh))
            raise ValueError

        # check if input mass is in mass nodes
        m = np.abs(self._mass_nodes - mass)<1e-3
        if m.sum()>0:
            alm = (self._alm1, self._alm2)[mass<=1.1]
            # get track
            track = self._load_track(y, feh, mass, alm)
            # put track into cache
            trackid = self._get_trackid(y, feh, mass, alm)
            self._track_data[trackid] = track

            if minage!=0:
                m = track[2] > minage
                track = (track[0][m], track[1][m], track[2][m],
                         track[3][m], track[4][m])
            if n>0:
                track = interpolate_data(track, n)
        else:
            track_lst = []
            imass = self._get_inodes(self._mass_nodes, mass)
            mass_lst = self._mass_nodes[imass:imass+4]
            for _mass in mass_lst:
                alm = (self._alm1, self._alm2)[_mass<=1.1]
                trackid = self._get_trackid(y, feh, _mass, alm)
                track = self._track_data[trackid]
                m = track[2] > minage
                track = (track[0][m], track[1][m], track[2][m],
                         track[3][m], track[4][m])
                track = interpolate_data(track, n)
                track_lst.append(track)
            track = interpolate_param(track_lst, mass_lst, mass)
        return track

    def _get_track_of_feh(self, y, feh, mass, n=0, minage=0):
        if y not in self._y_nodes:
            print('Error: Y = %g not in y_nodes')
            raise ValueError

        if feh in self._feh_nodes:
            # load missing tracks
            #trackid0 = self._get_trackid(y, feh, 1.0, self._alm1)
            #if trackid0 not in self._track_data:
            #    self._load_track(y, feh, mass, alm)
            track = self._get_track_of_mass(y, feh, mass, n, minage)
        else:
            # need interpolation over feh
            track_lst = []
            ifeh = self._get_inodes(self._feh_nodes, feh)
            feh_lst = self._feh_nodes[ifeh:ifeh+4]
            for _feh in feh_lst:
                # load missing tracks
                #trackid0 = self._get_trackid(y, _feh, 1.0, self._alm1)
                #if trackid0 not in self._track_data:
                #    self._load_track(y, _feh)
                track = self._get_track_of_mass(y, _feh, mass, n, minage)
                track_lst.append(track)
            track = interpolate_param(track_lst, feh_lst, feh)
        return track

YaPSI = _YaPSI()
