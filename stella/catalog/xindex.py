#!/usr/bin/env python
import os
import numpy as np

from ..utils.asciifile import find_sortedfile, quickfind_sortedfile

xindex_path = os.path.join(os.getenv('STELLA_DATA'),'catalog/xindex')

def HIP_to_HD(starname):
    if starname[0:3]=='HIP':
        hip = int(starname[3:])

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

def HIP_to_BD(starname):
    if starname[0:3]=='HIP':
        hip = int(starname[3:])

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

def HIP_to_CD(starname):
    if starname[0:3]=='HIP':
        hip = int(starname[3:])

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

def HIP_to_TYC(starname):
    if starname[0:3]=='HIP':
        hip = int(starname[3:])

    fn = '%s/HIP-TYC.fits'%xindex_path
    f = pf.open(fn)
    data = f[1].data
    f.close()
    m = data['HIP']==hip
    if m.sum()==0:
        return None
    else:
        return ['TYC %d-%d-%d'%(rec['TYC1'],rec['TYC2'],rec['TYC3']) for rec in data[m]]

def HIP_to_2MASS(starname,full=False):
    if starname[0:3]=='HIP':
        hip = int(starname[3:])

    fn = '%s/HIP-2MASS.fits'%xindex_path
    f = pf.open(fn)
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

def HD_to_HIP(starname):
    from .name import get_regular_name
    starname = get_regular_name(starname)
    hd = starname[2:].strip()

    fn = '%s/HD-HIP.csv'%xindex_path
    f1 = lambda row: row.split(',')[0].strip()
    f2 = lambda row: 'HIP %d'%(int(row.split(',')[1]))
    HIPname = quickfind_sortedfile(hd,fn,f1,f2)

    if HIPname == None:
        return None
    else:
        return [HIPname]

def HD_to_TYC(starname):
    from .name import get_regular_name
    starname = get_regular_name(starname)
    hd = starname[2:].strip()

    fn = '%s/HD-TYC.csv'%xindex_path
    f1 = lambda row: row.split(',')[0].strip()
    f2 = lambda row: 'TYC %s'%(row.split(',')[1].strip())

    TYCname = quickfind_sortedfile(hd,fn,f1,f2)
    if TYCname == None:
        return None
    else:
        return [TYCname]

def BD_to_HIP(starname):
    from .name import get_regular_name
    starname = get_regular_name(starname)
    bd = starname[2:].strip()

    fn = '%s/BD-HIP.csv'%xindex_path
    f1 = lambda row: row.split(',')[0].strip()
    f2 = lambda row: 'HIP %d'%(int(row.split(',')[1]))

    HIPname = quickfind_sortedfile(bd,fn,f1,f2)
    if HIPname == None:
        return None
    else:
        return [HIPname]

def CD_to_HIP(starname):
    from .name import get_regular_name
    starname = get_regular_name(starname)
    cd = starname[2:].strip()

    fn = '%s/CD-HIP.csv'%xindex_path
    f1 = lambda row: row.split(',')[0].strip()
    f2 = lambda row: 'HIP %d'%(int(row.split(',')[1]))
    
    HIPname = quickfind_sortedfile(cd,fn,f1,f2)
    if HIPname == None:
        return None
    else:
        return [HIPname]

def TYC_to_HIP(starname):
    if starname[0:3]=='TYC':
        g = starname[3:].split('-')
        tyc1,tyc2,tyc3 = int(g[0]),int(g[1]),int(g[2])

    fn = '%s/TYC-HIP.fits'%xindex_path
    f = pf.open(fn)
    data = f[1].data
    f.close()
    m1 = data['TYC1']==tyc1
    m2 = data[m1]['TYC2']==tyc2
    m3 = data[m1][m2]['TYC3']==tyc3
    if m3.sum()==0:
        return None
    else:
        return ['HIP %d'%rec['HIP'] for rec in data[m1][m2][m3]]

def TYC_to_2MASS(starname,full=False):
    if starname[0:3]=='TYC':
        g = starname[3:].split('-')
        tyc1,tyc2,tyc3 = int(g[0]),int(g[1]),int(g[2])

    fn = '%s/TYC-2MASS.fits'%xindex_path
    f = pf.open(fn)
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

def G_to_TYC(starname):
    if starname[0:2]=='G ':
        Gname = starname[2:].strip()

    fn = '%s/G-TYC.csv'%xindex_path
    f1 = lambda row: row.split(',')[0].strip()
    f2 = lambda row: 'TYC %s'%(row.split(',')[1].strip())

    TYCname = quickfind_sortedfile(Gname,fn,f1,f2)
    if TYCname == None:
        return None
    else:
        return [TYCname]

def KIC_to_KOI(string):
    fn = '%s/KIC-KOI.csv'%xindex_path
    kic = int(string)
    f1 = lambda row: int(row.split(',')[0])
    f2 = lambda row: int(row.split(',')[1])
    return quickfind_sortedfile(kic,fn,f1,f2)

def KIC_to_Kepler(string):
    fn = '%s/KIC-Kepler.csv'%xindex_path
    kic = int(string)
    f1 = lambda row: int(row.split(',')[0])
    f2 = lambda row: int(row.split(',')[1])
    return quickfind_sortedfile(kic,fn,f1,f2)

def KOI_to_KIC(string):
    fn = '%s/KOI-KIC.csv'%xindex_path
    koi = int(string)
    f1 = lambda row: int(row.split(',')[0])
    f2 = lambda row: int(row.split(',')[1])
    if koi < 100:
        return find_sortedfile(koi,fn,f1,f2)
    else:
        return quickfind_sortedfile(koi,fn,f1,f2)

def KOI_to_Kepler(string):
    fn = '%s/KOI-Kepler.csv'%xindex_path
    koi = int(string)
    f1 = lambda row: int(row.split(',')[0])
    f2 = lambda row: int(row.split(',')[1])
    if koi < 100:
        return find_sortedfile(koi,fn,f1,f2)
    else:
        return quickfind_sortedfile(koi,fn,f1,f2)

def Kepler_to_KIC(string):
    fn = '%s/Kepler-KIC.csv'%xindex_path
    kepler = int(string)
    f1 = lambda row: int(row.split(',')[0])
    f2 = lambda row: int(row.split(',')[1])
    if kepler < 100:
        return find_sortedfile(kepler,fn,f1,f2)
    else:
        return quickfind_sortedfile(kepler,fn,f1,f2)

def Kepler_to_KOI(string):
    fn = '%s/Kepler-KOI.csv'%xindex_path
    kepler = int(string)
    f1 = lambda row: int(row.split(',')[0])
    f2 = lambda row: int(row.split(',')[1])
    if kepler < 100:
        return find_sortedfile(kepler,fn,f1,f2)
    else:
        return quickfind_sortedfile(kepler,fn,f1,f2)

