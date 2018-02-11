#!/usr/bin/env python
import re
import numpy as np
from . import xindex

def get_catalog(starname):
    '''
    Return name of the star catalog from name of the star.
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

def cross_starnames(starname):
    name_lst = {}
    cat = get_catalog(starname)
    if cat == None:
        return None

    if cat == 'HIP':
        name_lst['HIP'] = [_get_regular_HIP_name(starname)]
        name_lst['HD'] = xindex.HIP_to_HD(name_lst['HIP'][0])
        name_lst['BD'] = xindex.HIP_to_BD(name_lst['HIP'][0])
        name_lst['CD'] = xindex.HIP_to_CD(name_lst['HIP'][0])
        name_lst['TYC'] = xindex.HIP_to_TYC(name_lst['HIP'][0])

        # fix two TYC for one HIP
        if name_lst['TYC']!=None and len(name_lst['TYC'])>1 and name_lst['HD']!=None:
            tmp = xindex.HD_to_TYC(name_lst['HD'][0])
            if len(tmp)==1 and tmp[0] in name_lst['TYC']:
                name_lst['TYC'] = tmp

        # if HIP->TYC failed, find TYC by HD
        if name_lst['TYC'] == None and name_lst['HD']!=None:
            name_lst['TYC'] = xindex.HD_to_TYC(name_lst['HD'][0])

        name_lst['2MASS'] = xindex.HIP_to_2MASS(name_lst['HIP'][0])
    elif cat == 'HD':
        name_lst['HD'] = [_get_regular_HD_name(starname)]
        name_lst['HIP'] = xindex.HD_to_HIP(name_lst['HD'][0])
        if name_lst['HIP'] == None:
            # no HIP name
            name_lst['TYC'] = xindex.HD_to_TYC(name_lst['HD'][0])
            if name_lst['TYC'] != None:
                name_lst['2MASS'] = xindex.TYC_to_2MASS(name_lst['TYC'][0])
        else:
            # if has HIP name
            name_lst['BD'] = xindex.HIP_to_BD(name_lst['HIP'][0])
            name_lst['CD'] = xindex.HIP_to_CD(name_lst['HIP'][0])
            name_lst['TYC'] = xindex.HIP_to_TYC(name_lst['HIP'][0])

            # fix two TYC for one HIP
            if name_lst['TYC']!=None and len(name_lst['TYC'])>1 and name_lst['HD']!=None:
                tmp = xindex.HD_to_TYC(name_lst['HD'][0])
                if len(tmp)==1 and tmp[0] in name_lst['TYC']:
                    name_lst['TYC'] = tmp

            # if HIP->TYC failed, find TYC by HD
            if name_lst['TYC'] == None:
                name_lst['TYC'] = xindex.HD_to_TYC(name_lst['HD'][0])

            name_lst['2MASS'] = xindex.HIP_to_2MASS(name_lst['HIP'][0])

    elif cat == 'BD':
        name_lst['BD'] = [_get_regular_BD_name(starname)]
        name_lst['HIP'] = xindex.BD_to_HIP(name_lst['BD'][0])
        if name_lst['HIP'] != None:
            name_lst['HD'] = xindex.HIP_to_HD(name_lst['HIP'][0])
            name_lst['CD'] = xindex.HIP_to_CD(name_lst['HIP'][0])
            name_lst['TYC'] = xindex.HIP_to_TYC(name_lst['HIP'][0])

            # fix two TYC for one HIP
            if name_lst['TYC']!=None and len(name_lst['TYC'])>1 and name_lst['HD']!=None:
                tmp = xindex.HD_to_TYC(name_lst['HD'][0])
                if len(tmp)==1 and tmp[0] in name_lst['TYC']:
                    name_lst['TYC'] = tmp

            name_lst['2MASS'] = xindex.HIP_to_2MASS(name_lst['HIP'][0])
    elif cat == 'CD':
        name_lst['CD'] = [_get_regular_CD_name(starname)]
        name_lst['HIP'] = xindex.CD_to_HIP(name_lst['CD'][0])
        if name_lst['HIP'] != None:
            name_lst['HD']    = xindex.HIP_to_HD(name_lst['HIP'][0])
            name_lst['BD']    = xindex.HIP_to_BD(name_lst['HIP'][0])
            name_lst['TYC']   = xindex.HIP_to_TYC(name_lst['HIP'][0])

            # fix two TYC for one HIP
            if name_lst['TYC']!=None and len(name_lst['TYC'])>1 and name_lst['HD']!=None:
                tmp = xindex.HD_to_TYC(name_lst['HD'][0])
                if len(tmp)==1 and tmp[0] in name_lst['TYC']:
                    name_lst['TYC'] = tmp

            name_lst['2MASS'] = xindex.HIP_to_2MASS(name_lst['HIP'][0])
    elif cat == 'G':
        name_lst['G'] = [_get_regular_G_name(starname)]
        name_lst['TYC'] = xindex.G_to_TYC(name_lst['G'][0])
        if name_lst['TYC'] != None:
            name_lst['2MASS'] = xindex.TYC_to_2MASS(name_lst['TYC'][0])
            name_lst['HIP']   = xindex.TYC_to_HIP(name_lst['TYC'][0])
            if name_lst['HIP'] != None:
                name_lst['HD'] = xindex.HIP_to_HD(name_lst['HIP'][0])
                name_lst['BD'] = xindex.HIP_to_BD(name_lst['HIP'][0])
                name_lst['CD'] = xindex.HIP_to_CD(name_lst['HIP'][0])
    elif cat == 'TYC':
        name_lst['TYC'] = [_get_regular_TYC_name(starname)]
        name_lst['2MASS'] = xindex.TYC_to_2MASS(name_lst['TYC'][0])
        name_lst['HIP']   = xindex.TYC_to_HIP(name_lst['TYC'][0])


    # delete those None names
    res_lst = {}
    for cat in name_lst:
        name = name_lst[cat]
        if name != None:
            res_lst[cat] = name
    return res_lst



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

