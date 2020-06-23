import os
import numpy as np
import astropy.io.fits as fits

from .name import get_regular_name, _get_HIP_number, _get_KIC_number
from ..utils.asciifile import find_sortedfile, quickfind_sortedfile

xindex_path = os.path.join(os.getenv('STELLA_DATA'), 'catalog/xindex')

def cross_starnames(starname):
    name_lst = {}
    cat = get_catalog(starname)
    if cat == None:
        return None

    if cat == 'HIP':
        name_lst['HIP'] = [_get_regular_HIP_name(starname)]
        name_lst['HD']  = HIP_to_HD(name_lst['HIP'][0])
        name_lst['BD']  = HIP_to_BD(name_lst['HIP'][0])
        name_lst['CD']  = HIP_to_CD(name_lst['HIP'][0])
        name_lst['TYC'] = HIP_to_TYC(name_lst['HIP'][0])

        # fix two TYC for one HIP
        if name_lst['TYC']!=None and len(name_lst['TYC'])>1 and name_lst['HD']!=None:
            tmp = HD_to_TYC(name_lst['HD'][0])
            if len(tmp)==1 and tmp[0] in name_lst['TYC']:
                name_lst['TYC'] = tmp

        # if HIP->TYC failed, find TYC by HD
        if name_lst['TYC'] == None and name_lst['HD']!=None:
            name_lst['TYC'] = HD_to_TYC(name_lst['HD'][0])

        name_lst['2MASS'] = HIP_to_2MASS(name_lst['HIP'][0])
    elif cat == 'HD':
        name_lst['HD'] = [_get_regular_HD_name(starname)]
        name_lst['HIP'] = HD_to_HIP(name_lst['HD'][0])
        if name_lst['HIP'] == None:
            # no HIP name
            name_lst['TYC'] = HD_to_TYC(name_lst['HD'][0])
            if name_lst['TYC'] != None:
                name_lst['2MASS'] = TYC_to_2MASS(name_lst['TYC'][0])
        else:
            # if has HIP name
            name_lst['BD']  = HIP_to_BD(name_lst['HIP'][0])
            name_lst['CD']  = HIP_to_CD(name_lst['HIP'][0])
            name_lst['TYC'] = HIP_to_TYC(name_lst['HIP'][0])

            # fix two TYC for one HIP
            if name_lst['TYC']!=None and len(name_lst['TYC'])>1 and name_lst['HD']!=None:
                tmp = HD_to_TYC(name_lst['HD'][0])
                if len(tmp)==1 and tmp[0] in name_lst['TYC']:
                    name_lst['TYC'] = tmp

            # if HIP->TYC failed, find TYC by HD
            if name_lst['TYC'] == None:
                name_lst['TYC'] = HD_to_TYC(name_lst['HD'][0])

            name_lst['2MASS'] = HIP_to_2MASS(name_lst['HIP'][0])

    elif cat == 'BD':
        name_lst['BD'] = [_get_regular_BD_name(starname)]
        name_lst['HIP'] = BD_to_HIP(name_lst['BD'][0])
        if name_lst['HIP'] != None:
            name_lst['HD'] = HIP_to_HD(name_lst['HIP'][0])
            name_lst['CD'] = HIP_to_CD(name_lst['HIP'][0])
            name_lst['TYC'] = HIP_to_TYC(name_lst['HIP'][0])

            # fix two TYC for one HIP
            if name_lst['TYC']!=None and len(name_lst['TYC'])>1 and name_lst['HD']!=None:
                tmp = HD_to_TYC(name_lst['HD'][0])
                if len(tmp)==1 and tmp[0] in name_lst['TYC']:
                    name_lst['TYC'] = tmp

            name_lst['2MASS'] = HIP_to_2MASS(name_lst['HIP'][0])
    elif cat == 'CD':
        name_lst['CD'] = [_get_regular_CD_name(starname)]
        name_lst['HIP'] = CD_to_HIP(name_lst['CD'][0])
        if name_lst['HIP'] != None:
            name_lst['HD']    = HIP_to_HD(name_lst['HIP'][0])
            name_lst['BD']    = HIP_to_BD(name_lst['HIP'][0])
            name_lst['TYC']   = HIP_to_TYC(name_lst['HIP'][0])

            # fix two TYC for one HIP
            if name_lst['TYC']!=None and len(name_lst['TYC'])>1 and name_lst['HD']!=None:
                tmp = HD_to_TYC(name_lst['HD'][0])
                if len(tmp)==1 and tmp[0] in name_lst['TYC']:
                    name_lst['TYC'] = tmp

            name_lst['2MASS'] = HIP_to_2MASS(name_lst['HIP'][0])
    elif cat == 'G':
        name_lst['G'] = [_get_regular_G_name(starname)]
        name_lst['TYC'] = G_to_TYC(name_lst['G'][0])
        if name_lst['TYC'] != None:
            name_lst['2MASS'] = TYC_to_2MASS(name_lst['TYC'][0])
            name_lst['HIP']   = TYC_to_HIP(name_lst['TYC'][0])
            if name_lst['HIP'] != None:
                name_lst['HD'] = HIP_to_HD(name_lst['HIP'][0])
                name_lst['BD'] = HIP_to_BD(name_lst['HIP'][0])
                name_lst['CD'] = HIP_to_CD(name_lst['HIP'][0])
    elif cat == 'TYC':
        name_lst['TYC'] = [_get_regular_TYC_name(starname)]
        name_lst['2MASS'] = TYC_to_2MASS(name_lst['TYC'][0])
        name_lst['HIP']   = TYC_to_HIP(name_lst['TYC'][0])


    # delete those None names
    res_lst = {}
    for cat in name_lst:
        name = name_lst[cat]
        if name != None:
            res_lst[cat] = name
    return res_lst


def HIP_to_HD(name):
    """Convert an HIP name in *Hipparcos Catalogue* to HD name in *Henry Draper
    Catalogue*.

    Args:
        name (str or int): Name of star in *Hipparcos Catalogue*.
    """

    hip = _get_HIP_number(name)

    filename = os.path.join(xindex_path, 'HIP-HD.csv')
    f1 = lambda row: int(row.split(',')[0])
    f2 = lambda row: 'HD '+row.split(',')[1].strip()
    if hip<100:
        HDname = find_sortedfile(hip, filename, f1, f2)
    else:
        HDname = quickfind_sortedfile(hip, filename, f1, f2)

    if HDname == None:
        return None
    else:
        return [HDname]

def HIP_to_BD(name):
    """Convert an HIP name in *Hipparcos Catalogue* to BD name in *Bonner
    Durchmusterung*.

    Args:
        name (str or int): Name of star in *Hipparcos Catalogue*.

    See also:
        * :ref:`catalog_hip`
        * :func:`stella.catalog.xindex.HIP_to_2MASS`
        * :func:`stella.catalog.xindex.HIP_to_CD`
        * :func:`stella.catalog.xindex.HIP_to_Gaia`
        * :func:`stella.catalog.xindex.HIP_to_HD`
        * :func:`stella.catalog.xindex.HIP_to_TYC`
    """

    hip = _get_HIP_number(name)

    filename = os.path.join(xindex_path, 'HIP-BD.csv')
    f1 = lambda row: int(row.split(',')[0])
    f2 = lambda row: 'BD '+row.split(',')[1].strip()
    if hip<100:
        HIPname = find_sortedfile(hip, filename, f1, f2)
    else:
        HIPname = quickfind_sortedfile(hip, filename, f1, f2)

    if HIPname == None:
        return None
    else:
        return [HIPname]

def HIP_to_CD(name):
    """Convert an HIP name in *Hipparcos Catalogue* to CD name in *Cordoba
    Durchmusterung*.

    Args:
        name (str or int): Name of star in *Hipparcos Catalogue*.

    See also:
        * :ref:`catalog_hip`
        * :func:`stella.catalog.xindex.HIP_to_HD`
        * :func:`stella.catalog.xindex.HIP_to_2MASS`
        * :func:`stella.catalog.xindex.HIP_to_BD`
        * :func:`stella.catalog.xindex.HIP_to_Gaia`
        * :func:`stella.catalog.xindex.HIP_to_HD`
        * :func:`stella.catalog.xindex.HIP_to_TYC`
    """

    hip = _get_HIP_number(name)

    filename = os.path.join(xindex_path, 'HIP-CD.csv')
    f1 = lambda row: int(row.split(',')[0])
    f2 = lambda row: 'CD '+row.split(',')[1].strip()
    if hip<100:
        CDname = find_sortedfile(hip, filename, f1, f2)
    else:
        CDname = quickfind_sortedfile(hip, filename, f1, f2)

    if CDname == None:
        return None
    else:
        return [CDname]

def HIP_to_TYC(name):
    """Convert an HIP name in *Hipparcos Catalogue* to TYC name in *Tycho-2
    Catalogue*.

    Args:
        name (str or int): Name of star in *Hipparcos Catalogue*.

    See also:
        * :ref:`catalog_hip`
        * :func:`stella.catalog.xindex.HIP_to_2MASS`
        * :func:`stella.catalog.xindex.HIP_to_BD`
        * :func:`stella.catalog.xindex.HIP_to_CD`
        * :func:`stella.catalog.xindex.HIP_to_Gaia`
        * :func:`stella.catalog.xindex.HIP_to_HD`
    """

    hip = _get_HIP_number(name)

    filename = os.path.join(xindex_path, 'HIP-TYC.fits')
    f = fits.open(filename)
    data = f[1].data
    f.close()
    m = data['HIP']==hip
    if m.sum()==0:
        return None
    else:
        return ['TYC %d-%d-%d'%(rec['TYC1'],rec['TYC2'],rec['TYC3']) for rec in data[m]]

def HIP_to_2MASS(name, full=False):
    """Convert an HIP name in *Hipparcos Catalogue* to 2MASS name.

    Args:
        name (str or int): Name of star in *Hipparcos Catalogue*.
        full (bool):

    See also:
        * :ref:`catalog_hip`
        * :func:`stella.catalog.xindex.HIP_to_CD`
        * :func:`stella.catalog.xindex.HIP_to_Gaia`
        * :func:`stella.catalog.xindex.HIP_to_HD`
        * :func:`stella.catalog.xindex.HIP_to_TYC`
    """

    hip = _get_HIP_number(name)

    filename = os.path.join(xindex_path, 'HIP-2MASS.fits')
    f = fits.open(filename)
    data = f[1].data
    i = np.searchsorted(data['HIP'],hip)
    row = data[i]
    des_2mass = '2MASS J%s'%row['2MASS']
    f.close()
    if full:
        return [(des_2mass, row['Jmag'], row['Hmag'], row['Kmag'],
                row['e_Jmag'], row['e_Hmag'], row['e_Kmag'])]
    else:
        return [des_2mass]

def HD_to_HIP(name):
    """Convert an HD name in *Henry Draper Catalogue* to HIP name in *Hipparcos
    Catalogue*.

    Args:
        name (str or int): Name of star in *Henry Draper Catalogue*.
    """
    name = get_regular_name(name)
    hd = name[2:].strip()

    filename = os.path.join(xindex_path, 'HD-HIP.csv')
    f1 = lambda row: row.split(',')[0].strip()
    f2 = lambda row: 'HIP %d'%(int(row.split(',')[1]))
    HIPname = quickfind_sortedfile(hd, filename, f1, f2)

    if HIPname is None:
        return None
    else:
        return [HIPname]

def HD_to_TYC(name):
    """Convert an HD name in *Henry Draper Catalogue* to TYC name in *Tycho-2
    Catalogue*.

    Args:
        name (str or int): Name of star in *Henry Draper Catalogue*.
    """
    name = get_regular_name(name)
    hd = name[2:].strip()

    filename = os.path.join(xindex_path, 'HD-TYC.csv')
    f1 = lambda row: row.split(',')[0].strip()
    f2 = lambda row: 'TYC %s'%(row.split(',')[1].strip())

    TYCname = quickfind_sortedfile(hd, filename, f1, f2)
    if TYCname == None:
        return None
    else:
        return [TYCname]

def BD_to_HIP(name):
    """Convert a BD name in *Bonner Durchmusterung* to HIP name in *Hipparcos
    Catalogue*.

    Args:
        name (str): Name of star in *Bonner Durchmusterung*.
    """
    name = get_regular_name(name)
    bd = name[2:].strip()

    filename = os.path.join(xindex_path, 'BD-HIP.csv')
    f1 = lambda row: row.split(',')[0].strip()
    f2 = lambda row: 'HIP %d'%(int(row.split(',')[1]))

    HIPname = quickfind_sortedfile(bd, filename, f1, f2)
    if HIPname == None:
        return None
    else:
        return [HIPname]

def CD_to_HIP(name):
    """Convert a CD name in *Cordoba Durchmusterung* to HIP name in *Hipparcos
    Catalogue*.

    Args:
        name (str): Name of star in *Cordoba Durchmusterung*.
    """
    name = get_regular_name(name)
    cd = name[2:].strip()

    filename = os.path.join(xindex_path, 'CD-HIP.csv')
    f1 = lambda row: row.split(',')[0].strip()
    f2 = lambda row: 'HIP %d'%(int(row.split(',')[1]))
    
    HIPname = quickfind_sortedfile(cd, filename, f1, f2)
    if HIPname == None:
        return None
    else:
        return [HIPname]

def TYC_to_HIP(name):
    """Convert a TYC name in *Tycho-2 Catalogue* to HIP name in *Hipparcos
    Catalogue*.
    """
    if name[0:3] == 'TYC':
        g = name[3:].split('-')
        tyc1, tyc2, tyc3 = int(g[0]), int(g[1]), int(g[2])

    filename = os.path.join(xindex_path, 'TYC-HIP.fits')
    f = fits.open(filename)
    data = f[1].data
    f.close()
    m1 = data['TYC1']==tyc1
    m2 = data[m1]['TYC2']==tyc2
    m3 = data[m1][m2]['TYC3']==tyc3
    if m3.sum()==0:
        return None
    else:
        return ['HIP %d'%rec['HIP'] for rec in data[m1][m2][m3]]

def TYC_to_2MASS(name, full=False):
    """Convert a TYC name to 2MASS name.

    Args:
        name (str):
        full (bool):

    Returns:

    """
    if name[0:3] == 'TYC':
        g = name[3:].split('-')
        tyc1, tyc2, tyc3 = int(g[0]), int(g[1]), int(g[2])

    filename = os.path.join(xindex_path, 'TYC-2MASS.fits')
    f = fits.open(filename)
    data = f[1].data
    f.close()

    m1 = data['TYC1']==tyc1
    m2 = data[m1]['TYC2']==tyc2
    m3 = data[m1][m2]['TYC3']==tyc3
    if m3.sum()==0:
        return None
    elif full:
        return [('2MASS J%s'%row['2MASS'], row['Jmag'], row['Hmag'],
                 row['Kmag'],row['e_Jmag'], row['e_Hmag'], row['e_Kmag']) for
                 row in data[m1][m2][m3]
                ]
    else:
        return ['2MASS J%s'%row['2MASS'] for row in data[m1][m2][m3]]

def G_to_TYC(name):
    """Convert a G name to TYC name in *Tycho-2 Catalogue*.

    Args:
        name (str):

    Returns:
        list
    """
    if name[0:2] == 'G ':
        Gname = name[2:].strip()

    filename = os.path.join(xindex_path, 'G-TYC.csv')
    f1 = lambda row: row.split(',')[0].strip()
    f2 = lambda row: 'TYC %s'%(row.split(',')[1].strip())

    TYCname = quickfind_sortedfile(Gname, filename, f1, f2)
    if TYCname == None:
        return None
    else:
        return [TYCname]

def KIC_to_KOI(name):
    """Convert a KIC name in *Kepler Input Catalog* to KOI name.

    Args:
        name (str):

    See also:
        * :func:`stella.catalog.xindex.KOI_to_KIC`
        * :func:`stella.catalog.xindex.KIC_to_Kepler`
        * :func:`stella.catalog.xindex.Kepler_to_KIC`
        * :func:`stella.catalog.xindex.Kepler_to_KOI`
        * :func:`stella.catalog.xindex.KOI_to_Kepler`
    """
    filename = os.path.join(xindex_path, 'KIC-KOI.csv')
    kic = _get_KIC_number(name)
    f1 = lambda row: int(row.split(',')[0])
    f2 = lambda row: int(row.split(',')[1])
    return quickfind_sortedfile(kic, filename, f1, f2)

def KIC_to_Kepler(name):
    """Convert a KIC name in *Kepler Input Catalog* to Kepler name.

    See also:
        * :func:`stella.catalog.xindex.Kepler_to_KIC`
        * :func:`stella.catalog.xindex.KIC_to_KOI`
        * :func:`stella.catalog.xindex.KOI_to_KIC`
        * :func:`stella.catalog.xindex.Kepler_to_KOI`
        * :func:`stella.catalog.xindex.KOI_to_Kepler`
    """
    filename = os.path.join(xindex_path, 'KIC-Kepler.csv')
    kic = _get_KIC_number(name)
    f1 = lambda row: int(row.split(',')[0])
    f2 = lambda row: int(row.split(',')[1])
    return quickfind_sortedfile(kic, filename, f1, f2)

def KOI_to_KIC(name):
    """Convert a KOI name to KIC name in *Kepler Input Catalog*.

    Args:
        name (str):

    See also:
        * :func:`stella.catalog.xindex.KIC_to_KOI`
        * :func:`stella.catalog.xindex.KOI_to_Kepler`
        * :func:`stella.catalog.xindex.Kepler_to_KOI`
        * :func:`stella.catalog.xindex.Kepler_to_KIC`
        * :func:`stella.catalog.xindex.KIC_to_Kepler`
    """
    filename = os.path.join(xindex_path, 'KOI-KIC.csv')
    koi = int(name)
    f1 = lambda row: int(row.split(',')[0])
    f2 = lambda row: int(row.split(',')[1])
    if koi < 100:
        return find_sortedfile(koi, filename, f1, f2)
    else:
        return quickfind_sortedfile(koi, filename, f1, f2)

def KOI_to_Kepler(name):
    """Convert a KOI name to Kepler name.

    Args:
        name (str):

    See also:
        * :func:`stella.catalog.xindex.Kepler_to_KOI`
        * :func:`stella.catalog.xindex.KOI_to_KIC`
        * :func:`stella.catalog.xindex.KIC_to_KOI`
        * :func:`stella.catalog.xindex.KIC_to_Kepler`
        * :func:`stella.catalog.xindex.Kepler_to_KIC`
    """
    filename = os.path.join(xindex_path, 'KOI-Kepler.csv')
    koi = int(name)
    f1 = lambda row: int(row.split(',')[0])
    f2 = lambda row: int(row.split(',')[1])
    if koi < 100:
        return find_sortedfile(koi, filename, f1, f2)
    else:
        return quickfind_sortedfile(koi, filename, f1, f2)

def Kepler_to_KIC(name):
    """Convert a Kepler name to KIC name.

    Args:
        name (str):

    See also:
        * :func:`stella.catalog.xindex.KIC_to_Kepler`
        * :func:`stella.catalog.xindex.Kepler_to_KOI`
        * :func:`stella.catalog.xindex.KOI_to_Kepler`
        * :func:`stella.catalog.xindex.KOI_to_KIC`
        * :func:`stella.catalog.xindex.KIC_to_KOI`
    """
    filename = os.path.join(xindex_path, 'Kepler-KIC.csv')
    kepler = int(name)
    f1 = lambda row: int(row.split(',')[0])
    f2 = lambda row: int(row.split(',')[1])
    if kepler < 100:
        return find_sortedfile(kepler, filename, f1, f2)
    else:
        return quickfind_sortedfile(kepler, filename, f1, f2)

def Kepler_to_KOI(name):
    """Convert a Kepler name to KOI name.

    Args:
        name (str):

    See also:
        * :func:`stella.catalog.xindex.KOI_to_Kepler`
        * :func:`stella.catalog.xindex.Kepler_to_KIC`
        * :func:`stella.catalog.xindex.KIC_to_Kepler`
        * :func:`stella.catalog.xindex.KOI_to_KIC`
        * :func:`stella.catalog.xindex.KIC_to_KOI`
    """
    filename = os.path.join(xindex_path, 'Kepler-KOI.csv')
    kepler = int(name)
    f1 = lambda row: int(row.split(',')[0])
    f2 = lambda row: int(row.split(',')[1])
    if kepler < 100:
        return find_sortedfile(kepler, filename, f1, f2)
    else:
        return quickfind_sortedfile(kepler, filename, f1, f2)

def HIP_to_Gaia(name):
    """Convert an HIP name in *Hipparcos Catalogue* to Gaia name.

    Args:
        name (str):
    """
    if isinstance(name, int):
        hip = name
    elif isinstance(name, str):
        if name[0:3]=='HIP':
            hip = int(name[3:])
        else:
            print('Unknown name: %s'%name)
    else:
        print('Unknown name')

    filename = os.path.join(xindex_path, 'HIP-Gaia.fits')
    data = fits.getdata()
    m = data['HIP']==hip
    if m.sum()==0:
        return None
    else:
        return data[m][0]['source_id']

def TOI_to_TIC(name):
    """Convert TOI (TESS Object of Interest) name to TIC name.

    Args:
        name (str):

    See also:
        * :func:`stella.catalog.xindex.TIC_to_TOI`
    """
    filename = os.path.join(xindex_path, 'TOI-TIC.csv')
    toi = int(name)
    f1 = lambda row: int(row.split(',')[0])
    f2 = lambda row: int(row.split(',')[1])
    return find_sortedfile(toi, filename, f1, f2)

def TIC_to_TOI(name):
    """Convert TIC name to TOI name.

    Args:
        name (str):

    See also:
        * :func:`stella.catalog.xindex.TOI_to_TIC`
    """
    filename = os.path.join(xindex_path, 'TIC-TOI.csv')
    tic = int(name)
    f1 = lambda row: int(row.split(',')[0])
    f2 = lambda row: int(row.split(',')[1])
    return find_sortedfile(tic, filename, f1, f2)
