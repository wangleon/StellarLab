import os
import numpy as np
import astropy.io.fits as fits

from .name import get_regular_name, _get_HIP_number, _get_KIC_number
from ..utils.asciifile import find_sortedfile, quickfind_sortedfile

xindex_path = os.path.join(os.getenv('STELLA_DATA'),'catalog/xindex')

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
    '''Convert HIP name to HD name.
    '''

    hip = _get_HIP_number(name)

    fn = '%s/HIP-HD.csv'%xindex_path
    f1 = lambda row: int(row.split(',')[0])
    f2 = lambda row: 'HD '+row.split(',')[1].strip()
    if hip<100:
        HDname = find_sortedfile(hip,fn,f1,f2)
    else:
        HDname = quickfind_sortedfile(hip,fn,f1,f2)

    if HDname == None:
        return None
    else:
        return [HDname]

def HIP_to_BD(name):
    '''Convert HIP name to BD name.

    See also:
        * :ref:`catalog_hip`
        * :func:`stella.catalog.xindex.HIP_to_2MASS`
        * :func:`stella.catalog.xindex.HIP_to_CD`
        * :func:`stella.catalog.xindex.HIP_to_Gaia`
        * :func:`stella.catalog.xindex.HIP_to_HD`
        * :func:`stella.catalog.xindex.HIP_to_TYC`
    '''

    hip = _get_HIP_number(name)

    fn = '%s/HIP-BD.csv'%xindex_path
    f1 = lambda row: int(row.split(',')[0])
    f2 = lambda row: 'BD '+row.split(',')[1].strip()
    if hip<100:
        HIPname = find_sortedfile(hip,fn,f1,f2)
    else:
        HIPname = quickfind_sortedfile(hip,fn,f1,f2)

    if HIPname == None:
        return None
    else:
        return [HIPname]

def HIP_to_CD(name):
    '''Convert HIP name to CD name.

    See also:
        * :ref:`catalog_hip`
        * :func:`stella.catalog.xindex.HIP_to_HD`
        * :func:`stella.catalog.xindex.HIP_to_2MASS`
        * :func:`stella.catalog.xindex.HIP_to_BD`
        * :func:`stella.catalog.xindex.HIP_to_Gaia`
        * :func:`stella.catalog.xindex.HIP_to_HD`
        * :func:`stella.catalog.xindex.HIP_to_TYC`
    '''

    hip = _get_HIP_number(name)

    fn = '%s/HIP-CD.csv'%xindex_path
    f1 = lambda row: int(row.split(',')[0])
    f2 = lambda row: 'CD '+row.split(',')[1].strip()
    if hip<100:
        CDname = find_sortedfile(hip,fn,f1,f2)
    else:
        CDname = quickfind_sortedfile(hip,fn,f1,f2)

    if CDname == None:
        return None
    else:
        return [CDname]

def HIP_to_TYC(name):
    '''Convert HIP name to TYC name.

    See also:
        * :ref:`catalog_hip`
        * :func:`stella.catalog.xindex.HIP_to_2MASS`
        * :func:`stella.catalog.xindex.HIP_to_BD`
        * :func:`stella.catalog.xindex.HIP_to_CD`
        * :func:`stella.catalog.xindex.HIP_to_Gaia`
        * :func:`stella.catalog.xindex.HIP_to_HD`
    '''

    hip = _get_HIP_number(name)

    fn = '%s/HIP-TYC.fits'%xindex_path
    f = fits.open(fn)
    data = f[1].data
    f.close()
    m = data['HIP']==hip
    if m.sum()==0:
        return None
    else:
        return ['TYC %d-%d-%d'%(rec['TYC1'],rec['TYC2'],rec['TYC3']) for rec in data[m]]

def HIP_to_2MASS(name,full=False):
    '''Convert HIP name to 2MASS name.

    See also:
        * :ref:`catalog_hip`
        * :func:`stella.catalog.xindex.HIP_to_CD`
        * :func:`stella.catalog.xindex.HIP_to_Gaia`
        * :func:`stella.catalog.xindex.HIP_to_HD`
        * :func:`stella.catalog.xindex.HIP_to_TYC`
    '''

    hip = _get_HIP_number(name)

    fn = '%s/HIP-2MASS.fits'%xindex_path
    f = fits.open(fn)
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
    '''Convert HD name to HIP name.
    '''
    name = get_regular_name(name)
    hd = name[2:].strip()

    fn = '%s/HD-HIP.csv'%xindex_path
    f1 = lambda row: row.split(',')[0].strip()
    f2 = lambda row: 'HIP %d'%(int(row.split(',')[1]))
    HIPname = quickfind_sortedfile(hd,fn,f1,f2)

    if HIPname == None:
        return None
    else:
        return [HIPname]

def HD_to_TYC(name):
    '''Convert HD name to TYC name.
    '''
    name = get_regular_name(name)
    hd = name[2:].strip()

    fn = '%s/HD-TYC.csv'%xindex_path
    f1 = lambda row: row.split(',')[0].strip()
    f2 = lambda row: 'TYC %s'%(row.split(',')[1].strip())

    TYCname = quickfind_sortedfile(hd,fn,f1,f2)
    if TYCname == None:
        return None
    else:
        return [TYCname]

def BD_to_HIP(name):
    '''Convert BD name to HIP name.
    '''
    name = get_regular_name(name)
    bd = name[2:].strip()

    fn = '%s/BD-HIP.csv'%xindex_path
    f1 = lambda row: row.split(',')[0].strip()
    f2 = lambda row: 'HIP %d'%(int(row.split(',')[1]))

    HIPname = quickfind_sortedfile(bd,fn,f1,f2)
    if HIPname == None:
        return None
    else:
        return [HIPname]

def CD_to_HIP(name):
    '''Convert CD name to HIP name.
    '''
    name = get_regular_name(name)
    cd = name[2:].strip()

    fn = '%s/CD-HIP.csv'%xindex_path
    f1 = lambda row: row.split(',')[0].strip()
    f2 = lambda row: 'HIP %d'%(int(row.split(',')[1]))
    
    HIPname = quickfind_sortedfile(cd,fn,f1,f2)
    if HIPname == None:
        return None
    else:
        return [HIPname]

def TYC_to_HIP(name):
    '''Convert TYC name to HIP name.
    '''
    if name[0:3]=='TYC':
        g = name[3:].split('-')
        tyc1,tyc2,tyc3 = int(g[0]),int(g[1]),int(g[2])

    fn = '%s/TYC-HIP.fits'%xindex_path
    f = fits.open(fn)
    data = f[1].data
    f.close()
    m1 = data['TYC1']==tyc1
    m2 = data[m1]['TYC2']==tyc2
    m3 = data[m1][m2]['TYC3']==tyc3
    if m3.sum()==0:
        return None
    else:
        return ['HIP %d'%rec['HIP'] for rec in data[m1][m2][m3]]

def TYC_to_2MASS(name,full=False):
    '''Convert TYC name to 2MASS name.
    '''
    if name[0:3]=='TYC':
        g = name[3:].split('-')
        tyc1,tyc2,tyc3 = int(g[0]),int(g[1]),int(g[2])

    fn = '%s/TYC-2MASS.fits'%xindex_path
    f = fits.open(fn)
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
    '''Convert G name to TYC name.
    '''
    if name[0:2]=='G ':
        Gname = name[2:].strip()

    fn = '%s/G-TYC.csv'%xindex_path
    f1 = lambda row: row.split(',')[0].strip()
    f2 = lambda row: 'TYC %s'%(row.split(',')[1].strip())

    TYCname = quickfind_sortedfile(Gname,fn,f1,f2)
    if TYCname == None:
        return None
    else:
        return [TYCname]

def KIC_to_KOI(name):
    '''Convert KIC name to KOI name.
    '''
    fn = '%s/KIC-KOI.csv'%xindex_path
    kic = _get_KIC_number(name)
    f1 = lambda row: int(row.split(',')[0])
    f2 = lambda row: int(row.split(',')[1])
    return quickfind_sortedfile(kic,fn,f1,f2)

def KIC_to_Kepler(name):
    '''Convert KIC name to Kepler name.
    '''
    fn = '%s/KIC-Kepler.csv'%xindex_path
    kic = _get_KIC_number(name)
    f1 = lambda row: int(row.split(',')[0])
    f2 = lambda row: int(row.split(',')[1])
    return quickfind_sortedfile(kic,fn,f1,f2)

def KOI_to_KIC(name):
    '''Convert KOI name to KIC name.
    '''
    fn = '%s/KOI-KIC.csv'%xindex_path
    koi = int(name)
    f1 = lambda row: int(row.split(',')[0])
    f2 = lambda row: int(row.split(',')[1])
    if koi < 100:
        return find_sortedfile(koi,fn,f1,f2)
    else:
        return quickfind_sortedfile(koi,fn,f1,f2)

def KOI_to_Kepler(name):
    '''Convert KOI name to Kepler name.
    '''
    fn = '%s/KOI-Kepler.csv'%xindex_path
    koi = int(name)
    f1 = lambda row: int(row.split(',')[0])
    f2 = lambda row: int(row.split(',')[1])
    if koi < 100:
        return find_sortedfile(koi,fn,f1,f2)
    else:
        return quickfind_sortedfile(koi,fn,f1,f2)

def Kepler_to_KIC(name):
    '''Convert Kepler name to KIC name.
    '''
    fn = '%s/Kepler-KIC.csv'%xindex_path
    kepler = int(name)
    f1 = lambda row: int(row.split(',')[0])
    f2 = lambda row: int(row.split(',')[1])
    if kepler < 100:
        return find_sortedfile(kepler,fn,f1,f2)
    else:
        return quickfind_sortedfile(kepler,fn,f1,f2)

def Kepler_to_KOI(name):
    '''
    Convert Kepler name to KOI name
    '''
    fn = '%s/Kepler-KOI.csv'%xindex_path
    kepler = int(name)
    f1 = lambda row: int(row.split(',')[0])
    f2 = lambda row: int(row.split(',')[1])
    if kepler < 100:
        return find_sortedfile(kepler,fn,f1,f2)
    else:
        return quickfind_sortedfile(kepler,fn,f1,f2)

def HIP_to_Gaia(name):
    '''Convert HIP name to Gaia name.
    '''
    if isinstance(name, int):
        hip = name
    elif isinstance(name, str):
        if name[0:3]=='HIP':
            hip = int(name[3:])
        else:
            print('Unknown name: %s'%name)
    else:
        print('Unknown name')

    data = fits.getdata('%s/HIP-Gaia.fits'%xindex_path)
    m = data['HIP']==hip
    if m.sum()==0:
        return None
    else:
        return data[m][0]['source_id']
         
