
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

def _get_TYC_number(starname):
    '''
    Private function to convert star name to TYC number (TYC1, TYC2, TYC3).
    '''
    if starname[0:3]=='TYC':
        g = starname[3:].split('-')
        tyc1, tyc2, tyc3 = int(g[0]), int(g[1]), int(g[2])
        return (tyc1, tyc2, tyc3)
    elif len(starname.split('-'))==3:
        g = starname.split('-')
        tyc1, tyc2, tyc3 = int(g[0]), int(g[1]), int(g[2])
        return (tyc1, tyc2, tyc3)
    else:
        return None

def _str_to_float(string):
    '''
    Convert string to float. Return *None* if failed.
    '''
    try:
        return float(string)
    except:
        return None

def _str_to_int(string):
    '''
    Convert string to integer. Return *None* if failed.
    '''
    try:
        return int(string)
    except:
        return None
