
def _get_star_number1(starname, key):
    '''
    Convert star name with the form of `SSS NNNN` to its integer number `NNNN`.

    Args:
        starname (integer or string): Name of a star
        key (string): Prefix of the star name
    Returns:
        integer: Number of the star in the catalog. If fail a *None* value will
            be returned.

    '''
    if isinstance(starname, int):
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
    '''Convert star name to an integer HIP number.

    Args:
        starname (string or integer): Name of the star
    Returns:
        integer: HIP number
    '''
    return _get_star_number1(starname, 'HIP')

def _get_KIC_number(starname):
    '''Convert star name to an integer KIC number.

    Args:
        starname (string or integer): Name of the star
    Returns:
        integer: KIC number
    '''
    return _get_star_number1(starname, 'KIC')

def _get_EPIC_number(starname):
    '''Convert star name to an integer EPIC number.

    Args:
        starname (string or integer): Name of the star
    Returns:
        integer: EPIC number
    '''
    return _get_star_number1(starname, 'EPIC')

def _get_TYC_number(starname):
    '''Convert star name to TYC number (TYC1, TYC2, TYC3).

    Args:
        starname (string): Name of the star
    Returns:
        tuple: A tuple of TYC numbers (TYC1, TYC2, TYC3)
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

def _str_to_float(string, exception_value=None):
    '''
    Convert string to float. Return `exception_value` if failed.
    '''
    try:
        return float(string)
    except:
        return exception_value

def _str_to_int(string, exception_value=None):
    '''
    Convert string to integer. Return `exception_value` if failed.
    '''
    try:
        return int(string)
    except:
        return exception_value
