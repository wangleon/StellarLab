import re
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table
from astroquery.simbad import Simbad

from . import name_parser

def get_names(names, catalogs=None):
    """Get names from Simbad data base.

    Args:
        names (list): A list of strings.
    """

    name_lst = {}

    for name in names:
        mobj = re.match('\*\s*([a-z]+\s+[a-zA-Z]{3})$', name)
        if mobj and name_parser.is_Bayer_name(mobj.group(1)):
            newname = name_parser.parse_Bayer_name(mobj.group(1))
            if 'Bayer' not in name_lst:
                name_lst['Bayer'] = []
            name_lst['Bayer'].append(newname)
            continue

        mobj = re.match('\*\s*(\d+\s+[a-zA-Z]{3})$', name)
        if mobj and name_parser.is_Flamsteed_name(mobj.group(1)):
            newname = name_parser.parse_Flamsteed_name(mobj.group(1))
            if 'Flamsteed' not in name_lst:
                name_lst['Flamsteed'] = []
            name_lst['Flamsteed'].append(newname)
            continue

        mobj = re.match('V\*\s*([\s\S]+)$', name)
        if mobj and name_parser.is_Var_name(mobj.group(1)):
            newname = name_parser.parse_Var_name(mobj.group(1))
            name_lst['Var'] = newname
            continue

        mobj = re.match('HD\s*(\d*[A-D]?)$', name)
        if mobj:
            name_lst['HD'] = mobj.group(1)
            continue

        mobj = re.match('HR\s*(\d+)$', name)
        if mobj:
            name_lst['HR'] = int(mobj.group(1))
            continue

        mobj = re.match('HIP\s*(\d+)$', name)
        if mobj:
            name_lst['HIP'] = int(mobj.group(1))
            continue

        if name_parser.is_DM_name(name):
            if 'DM' not in name_lst:
                name_lst['DM'] = []
            name_lst['DM'].append(name_parser.parse_DM_name(name))
            continue

        mobj = re.match('TYC\s*(\d*)\-(\d*)\-(\d)$', name)
        if mobj:
            name_lst['TYC'] = name_parser.parse_TYC_name(name)
            continue

        mobj = re.match('TIC\s*(\d*)$', name)
        if mobj:
            name_lst['TIC'] = int(mobj.group(1))
            continue

        mobj = re.match('GSC\s*(\d+\-\d+)$', name)
        if mobj:
            gsc = mobj.group(1)
            if len(gsc)>11:
                print('Warning: GSC name truncated:', gsc)
            name_lst['GSC'] = gsc
            continue

        mobj = re.match('CCDM\s*J(\S*)$', name)
        if mobj:
            if 'CCDM' in name_lst:
                print('CCDM already in names, ', name_lst['CCDM'], mobj.group(1))
            name_lst['CCDM'] = mobj.group(1)
            continue

        mobj = re.match('2MASS\s*J(\S*)$', name)
        if mobj:
            tmass = mobj.group(1)
            if len(tmass)>16:
                print('Warning: 2MASS name truncated:', tmass)
            name_lst['2MASS'] = tmass
            continue

        mobj = re.match('Gaia DR2 (\d*)', name)
        if mobj:
            gaia2 = int(mobj.group(1))
            name_lst['Gaia2'] = gaia2
            continue

        mobj = re.match('UCAC4\s*(\d+\-\d+)', name)
        if mobj:
            ucac4 = mobj.group(1).strip()
            if len(ucac4)>10:
                print('Warning: UCAC4 name truncated:', ucac4)
            name_lst['UCAC4'] = ucac4
            continue
    return name_lst

def select_main_name(name_lst):
    if 'Bayer' in name_lst:
        main_name = name_lst['Bayer']
        if isinstance(main_name, list):
            return main_name[0]
        else:
            return main_name
    elif 'Flamsteed' in name_lst:
        main_name = name_lst['Flamsteed']
        if isinstance(main_name, list):
            return main_name[0]
        else:
            return main_name
    elif 'Var' in name_lst:
        return name_lst['Var']
    elif 'HD' in name_lst:
        return 'HD {}'.format(name_lst['HD'])
    elif 'DM' in name_lst:
        main_name = name_lst['DM']
        if isinstance(main_name, list):
            return main_name[0]
        else:
            return main_name
    elif 'HIP' in name_lst:
        return 'HIP {}'.format(name_lst['HIP'])
    elif 'TYC' in name_lst:
        return 'TYC {}'.format(name_lst['TYC'])
    elif 'GSC' in name_lst:
        return 'GSC {}'.format(name_lst['GSC'])
    elif '2MASS' in name_lst:
        return '2MASS J{}'.format(name_lst['2MASS'])
    else:
        return ''

def make_crosstable(name_lst, coord_lst,
        columns=['HD','DM','HIP','TYC','TIC','GSC','CCDM','UCAC4','2MASS','Gaia2'],
        verbose = True,
        ):

    dtype = [
        ('name',        'S50'),
        ('RA_J2000',    'f8'),
        ('Dec_J2000',   'f8'),
        ]

    column_types = {
            'Bayer':    'S20',
            'Flamsteed':'S20',
            'Var':      'S20',
            'HD':       'S10',
            'DM':       'S20',
            'HR':       'i4',
            'HIP':      'i4',
            'TYC':      'S12',
            'TIC':      'i8',
            'GSC':      'S11',
            'UCAC4':    'S10',
            'CCDM':     'S15',
            '2MASS':    'S16',
            'Gaia2':    'i8',
            }
    for column in columns:
        if column in column_types:
            dtype.append((column, column_types[column]))

    crosstable = Table(dtype=dtype, masked=True)

    for name, coord in zip(name_lst, coord_lst):

        if verbose:
            print('searching around {:9.5f} {:9.5f} {}'.format(
                coord.ra.deg, coord.dec.deg, name))

        # query object ID
        result1 = Simbad.query_objectids(name)
        if result1 is None:
            print('\033[31mWarning: Missing object id:{}\033[0m'.format(name))
            continue
        crossname_lst = get_names(list(result1['ID']))

        # query object
        result2 = Simbad.query_object(name)
        if result2 is None:
            print('\033[31mWarning: Missing object:{}\033[0m'.format(name))
        row0 = result2[0]
        newcoord = SkyCoord(row0['RA'], row0['DEC'], unit=(u.hourangle, u.deg))

        ## query coordinate
        #for r in [3, 5, 10]:
        #    result3 = Simbad.query_region(coord, radius=r*u.arcsec)
        #    if result3 is not None:
        #        break
        #if result3 is None:
        #    print('\033[31mWarning: missing object in coord ({}, {}) {}\033[0m'.format(
        #        coord.ra.deg, coord.dec.deg, name))
        #    continue
        #row0 = result3[0]
        #newcoord = SkyCoord(row0['RA'], row0['DEC'], unit=(u.hourangle, u.deg))

        main_name = select_main_name(crossname_lst)

        item = [
            (main_name, False),
            (newcoord.ra.deg, False),
            (newcoord.dec.deg, False),
            ]

        for column in columns:
            if column in ['Bayer', 'Flamsteed']:
                if column in crossname_lst:
                    value = ';'.join(crossname_lst[column])
                    mask = False
                else:
                    value = ''
                    mask = True
            elif column in ['Var', 'HD', 'DM', 'TYC', 'GSC', 'CCDM',
                'UCAC4', '2MASS']:
                if column in crossname_lst:
                    value = crossname_lst[column]
                    if isinstance(value, list):
                        # if value is a list (e.g. DM), then multiple names
                        value = value[0]
                    mask = False
                else:
                    value = ''
                    mask = True
            elif column in ['HR', 'HIP', 'TIC', 'Gaia2']:
                if column in crossname_lst:
                    value = crossname_lst[column]
                    mask = False
                else:
                    value = 0
                    mask = True

            item.append((value, mask))

        value, mask = list(zip(*item))
        crosstable.add_row(value, mask=mask)


    crosstable['RA_J2000'].info.format  = '%9.5f'
    crosstable['Dec_J2000'].info.format = '%9.5f'
    for column in ['DM','CCDM']:
        if column in crosstable.colnames:
            lenlist = [len(row[column]) for row in crosstable
                        if row[column] is not np.ma.masked]
            maxlen = max(lenlist)
            crosstable[column].info.format = '<{}s'.format(maxlen)

    return crosstable
