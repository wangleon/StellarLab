
import numpy as np

def _get_star_number1(starname, key):
    '''
    Private function to convert star name to an integer HIP number.

    Parameters
    -----------
    starname : *integer*, *string*, *numpy.int_*
        Name of star
    key : *string*
        Prefix of star name

    Notes
    ------
    If fail, the `None` value will be returned.
    '''
    if isinstance(starname, int):
        return starname
    elif isinstance(starname, np.int_):
        return starname
    elif isinstance(starname, str):
        if starname[0:len(key)] == key:
            return int(starname[len(key):])
        elif starname.isdigit():
            return int(starname)
        else:
            return None
    else:
        return None

def _get_HIP_number(starname):
    '''
    Private function to convert star name to an integer HIP number.
    '''
    return _get_star_number1(starname, 'HIP')

def _get_KIC_number(starname):
    '''
    Private function to convert star name to an integer KIC number.
    '''
    return _get_star_number1(starname, 'KIC')
