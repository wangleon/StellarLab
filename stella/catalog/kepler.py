import os
from .base import _str_to_float, _str_to_int

planet_files = {
        1: 'ApJ.728.117.tablea1.dat', # Borucki et al. 2011a
        2: 'ApJ.736.19.table2.dat',   # Borucki et al. 2011b
        }


def load_systems(release):
    """Return a planetary system list in the given data releases.

    Args:
        release:
    """
    pass

def load_planets(release):
    """Return a planet list in the given data releases.

    Args:
        release:
    """
    pass

def find_system(koi, release):
    """Find parameters of a planetary system in the given data releases.

    Args:
        koi (int): KOI number of a planetary system.
        release (int): List of data releases.
    Return:
        dict: A dict containing parameters of planets in the system.
    """
    result = {}
    for dataset in release:
        if dataset == 1:
            pass
        elif dataset == 2:
            planet_file = planet_files[dataset]
            filename = os.path.join(os.getenv('STELLA_DATA'),
                        'catalog/Kepler/%s'%planet_file)
            infile = open(filename)
            for row in infile:
                record = _parse_planet_record_r2(row)
                if record['koi'] == koi:
                    if record['planet_id'] not in result:
                        result[record['planet_id']] = []
                    result[record['planet_id']].append(record)
                elif record['koi'] > koi:
                    break
            infile.close()
    return result

def find_planet(planet_id, release):
    """Find parameters of a planet in the given data releases.
    
    Args:
        planet_id (float): KOI number `NNN.NN` of a planet candidate.
        release (int): List of data releases.
    Return:
        list: A list containing parameter tuple as elements.
    """
    result = []
    for dataset in release:
        if dataset == 2:
            planet_file = planet_files[dataset]
            filename = os.path.join(os.getenv('STELLA_DATA'),
                        'catalog/Kepler/%s'%planet_file)
            infile = open(filename)
            for row in infile:
                record = _parse_planet_record_r2(row)
                if record['planet_id'] == planet_id:
                    result.append(record)
                    break
            infile.close()
    return result

def _parse_planet_record_r1(row):
    """Parse a planet record in the table of `Borucki+ 2011a
    <http://adsabs.harvard.edu/abs/2011ApJ...728..117B>`_.

    Args:
        row (str): The row in Borucki et al. 2011a.
    Returns:
        dict: A dict containing the planet parameters.
    """
    koi       = int(row[0:3])
    planet_id = float(row[0:6])
    r         = float(row[21:25])*11.209
    P         = float(row[34:41])
    rs        = float(row[53:58])
    rR        = r*6371/rs/695700.
    return {'koi': koi,  'planet_id': planet_id,
            'r'  : r,
            'P'  : P,
            'rR' : rR,
            }

def _parse_planet_record_r2(row):
    """Parse a planet record in the table of `Borucki+ 2011b
    <http://adsabs.harvard.edu/abs/2011ApJ...736...19B>`_.

    Args:
        row (str): The row in Borucki et al. 2011b.
    Returns:
        dict: A dict containing the planet parameters.
    """
    koi       = int(row[12:16])
    planet_id = float(row[12:19])
    Tdur      = _str_to_float(row[20:27])
    depth     = _str_to_float(row[28:34])
    P         = _str_to_float(row[60:73])
    P_err     = _str_to_float(row[74:86])
    aR        = _str_to_float(row[87:98])
    aR_err    = _str_to_float(row[99:110])
    rR        = _str_to_float(row[111:118])
    rR_err    = _str_to_float(row[119:126])
    b         = _str_to_float(row[127:133])
    b_err     = _str_to_float(row[134:139])
    r         = _str_to_float(row[140:145])
    a         = _str_to_float(row[146:151])
    Teq       = _str_to_int(row[152:156])

    return {
        'koi'   : koi,   'planet_id': planet_id,
        'Tdur'  : Tdur,
        'depth' : depth,
        'P'     : P,     'e_P'   : P_err,
        'r/R*'  : rR,    'e_r/R*': rR_err,
        'a/R*'  : aR,    'e_a/R*': aR_err,
        'b'     : b,     'e_b'   : b_err,
        'r'     : r,
        'a'     : a,
        'Teq'   : Teq,
        'status': 'candidate',
        'ref'   : 'Borucki et al. 2011b',
        }


def find_Kepler_cands_r2(koi):
    """Find planet candidates in the second release of *Kepler Mission* on Feb.
    2, 2011 (Borucki+ 2011).

    On Feb. 2, 2011, the *Kepler* team announced the second planet candidate
    list based on the data taken between May 2 and Sep. 16, 2010
    (`Borucki et al. 2011 <http://adsabs.harvard.edu/abs/2011ApJ...736...19B>`_).
    There are 1235 planet candidates associated with 997 host stars.

    .. csv-table:: Descriptions of returned parameters
        :header: Key, Type, Unit, Description
        :widths: 30, 30, 30, 100

        Tdur,   float32,   hour,        Transit duration
        depth,  float32,   ppm,         Transit depth
        P,      float32,   day,         Orbital period
        e_P,    float32,   day,         Uncertainty in orbital period
        r,      float32,   *R*:sub:`⊕`, Planet radius
        a/R*,   float32,   ,            Ratio of semi-major axis to stellar radius
        e_a/R*, float32,   ,            Uncertainty in a/R*
        r/R*,   float32,   ,            Ratio of planet radius to stellar radius
        e_r/R*, float32,   ,            Uncertainty in r/R*
        b,      float32,   ,            Impact parameter of transit
        e_b,    float32,   ,            Uncertainty in impact parameter
        a,      float32,   AU,          Semi-major axis
        Teq,    integer16, K,           Equilibrium temperature of planet
        status, string,    ,            Status. fixed to 'candidate'
        ref,    string,    ,            Reference. fixed to 'Borucki et al. 2011b'

    Args:
        koi (int or float): KOI number of the star (if finding all planets
            orbiting a star) or planet name (if finding a planet).
    Returns:
        dict: A dict containing parameters of all planets around a host star or
            the parameters of the input planet.

    """

    filename = os.path.join(os.getenv('STELLA_DATA'),
                'catalog/Kepler/ApJ.736.19.table2.dat')
    infile = open(filename)

    planet_lst = {}

    if isinstance(koi, int):
        # input is a KOI number of the host star
        obj = 's'
    elif isinstance(koi, float):
        # input is a planet id
        obj = 'p'
    else:
        # unrecognized name
        raise ValueError

    for row in infile:
        _koi_s = int(row[12:16])
        _koi_p = float(row[12:19])
        if (obj=='s' and _koi_s==koi) or (obj=='p' and _koi_p==koi):
            Tdur      = _str_to_float(row[20:27])
            depth     = _str_to_float(row[28:34])
            P         = _str_to_float(row[60:73])
            P_err     = _str_to_float(row[74:86])
            aR        = _str_to_float(row[87:98])
            aR_err    = _str_to_float(row[99:110])
            rR        = _str_to_float(row[111:118])
            rR_err    = _str_to_float(row[119:126])
            b         = _str_to_float(row[127:133])
            b_err     = _str_to_float(row[134:139])
            r         = _str_to_float(row[140:145])
            a         = _str_to_float(row[146:151])
            Teq       = _str_to_int(row[152:156])

            planet_data = {
                'Tdur'  : Tdur,
                'depth' : depth,
                'P'     : P,     'e_P'   : P_err,
                'r/R*'  : rR,    'e_r/R*': rR_err,
                'a/R*'  : aR,    'e_a/R*': aR_err,
                'b'     : b,     'e_b'   : b_err,
                'r'     : r,
                'a'     : a,
                'Teq'   : Teq,
                'status': candidate,
                'ref'   : 'Borucki et al. 2011b',
                }
            planet_lst[planet_id] = planet_data

        elif _koi_s > int(koi):
            break
        else:
            continue
    infile.close()

    if obj=='s':
        if len(planet_lst)>0:
            return planet_lst
        else:
            return None
    elif obj=='p':
        if len(planet_lst)==1:
            return planet_lst[planet_lst.keys()[0]]
        else:
            return None
    else:
        raise ValueError
