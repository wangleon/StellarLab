import re

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

def get_catalog(starname):
    '''Return the name of the star catalog from the name of star.
    '''
    starcat_re = {
        '^[Hh][Dd][\d\s]+[ABC]?$'          : 'HD',
        '^[Hh][Ii][Pp][\d\s]+$'            : 'HIP',
        '^[Hh][Rr][\d\s]+$'                : 'HR',
        '^[Bb][Dd][\+\-\d\s]+[a-zA-Z]?$'   : 'BD',
        '^[Cc][Dd][\-\d\s]+[a-zA-Z]?$'     : 'CD',
        '^TYC\s*\d+\-\d+\-\d?$'            : 'TYC',
        '^SAO[\d\s]+$'                     : 'SAO',
        '^GC[\d\s]+$'                      : 'GC',
        '^GCRV[\d\s]+$'                    : 'GCRV',
        '^G[\-\d\s]+$'                     : 'G',
        '^NLTT[\d\s]+$'                    : 'NLTT',
        '^LHS[\d\s]+[ABC]?$'               : 'LHS',
        '^LSPM[\d\s]+[NSEW]?$'             : 'LSPM',
        '^FK5[\d\s]+$'                     : 'FK5'
    }
    for exp in starcat_re:
        if re.match(exp, starname) != None:
            return starcat_re[exp]
    return None

def get_regular_name(starname):
    '''
    Get regular name of a star
    '''
    cat = get_catalog(starname)
    if cat == 'HIP':
        return _get_regular_HIP_name(starname)
    elif cat == 'HD':
        return _get_regular_HD_name(starname)
    elif cat == 'BD':
        return _get_regular_BD_name(starname)
    elif cat == 'CD':
        return _get_regular_CD_name(starname)
    elif cat == 'G':
        return _get_regular_G_name(starname)
    elif cat == 'TYC':
        return _get_regular_TYC_name(starname)
    else:
        return starname

def _get_regular_HIP_name(starname):
    '''Convert an HIP name to its regular form.
    '''

    if isinstance(starname, str):
        starname = starname.strip()
        starname = 'HIP '+starname[3:].strip()

    elif isinstance(starname, int):
        starname = 'HIP '+str(starname)

    else:
        raise ValueError

    return starname

def _get_regular_HD_name(starname):
    '''Convert an HD name to its regular form.
    '''

    if isinstance(starname, str):

        starname = starname.strip()

        if starname[-1].isalpha():
            comp = starname[-1]
            starname = 'HD %d %s'%(int(starname[2:-1]), comp)
        else:
            starname = 'HD %d'%(int(starname[2:]))

    elif isinstance(starname, int):
        starname = 'HD %d'%starname

    else:
        raise ValueError

    return starname


def _get_regular_BD_name(starname):
    '''Convert a BD name to its regular form.
    '''

    starname = starname.strip()

    # parse 'BD+33 23' to 'BD +33 23'
    starname = 'BD '+starname[2:].strip()

    # parse 'BD 33 23' to 'BD +33 23'
    g = starname.split()
    if g[1][0].isdigit():
        g[1] = '+'+g[1]
        starname = ' '.join(g)

    # parse 'BD +33 23A' to 'BD +33 23 A'
    # but keep 'BD +33 23a'
    if starname[-1].isalpha():
        comp = starname[-1]
        if comp.isupper():
            starname = starname[:-1].strip()+' '+comp
        elif comp.islower():
            starname = starname[:-1].strip()+comp
        else:
            raise ValueError

    # parse 'BD -3 23' to 'BD -03 23'
    g = starname.split()
    if len(g[1])!=3:
        pm = g[1][0]
        num = abs(int(g[1]))
        g[1] = pm + str(num).rjust(2,'0')

    # parse 'BD -03 0023' to 'BD -03 23'
    g[2] = str(int(g[2]))
    starname = ' '.join(g)

    return starname

def _get_regular_CD_name(starname):
    '''Convert a CD name to its regular form.
    '''

    starname = starname.strip()

    # parse 'CD-33 23' to 'CD -33 23'
    starname = 'CD '+starname[2:].strip()

    # parse 'CD -33 23A' to 'CD -33 23 A'
    if starname[-1].isalpha():
        comp = starname[-1]
        if comp.isupper():
            starname = starname[:-1].strip()+' '+comp
        elif comp.islower():
            starname = starname[:-1].strip()+comp
        else:
            raise ValueError

    # parse 'CD -22 0023' to 'CD -22 23'
    g = starname.split()
    g[2] = str(int(g[2]))
    starname = ' '.join(g)

    return starname


def _get_regular_G_name(starname):
    '''Convert a G name to its regular form.
    '''

    starname = starname.strip()
    if starname[0]=='G' and not starname[1].isalpha():
        starname = 'G '+starname[1:].strip()
    else:
        return None

    if starname[-1].isalpha():
        comp = starname[-1]
        starname = starname[0:-1].strip()+' '+comp
    else:
        comp = None

    # parse 'G 064-012' to 'G 64-12'
    g = starname.split()
    k = g[1].split('-')
    s = str(int(k[0]))+'-'+str(int(k[1]))
    if comp==None:
        starname = 'G '+s
    else:
        starname = 'G '+s+' '+comp

    return starname

def _get_regular_TYC_name(*args):
    '''Convert a TYC name to its regular form.
    '''
    if len(args)==1 and isinstance(args[0], str):
        starname = args[0].strip()
        starname = 'TYC '+starname[3:].strip()

    elif len(args)==3:
        starname = 'TYC '+'-'.join([str(args[0]),str(args[1]),str(args[2])])
    else:
        return None

    return starname

