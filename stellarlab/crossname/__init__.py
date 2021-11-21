import re
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as units
from astropy.table import Table
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier

from . import name_parser

def get_names(names, catalogs=None):
    """Get names from Simbad data base.

    Args:
        names (list): A list of strings.
    """

    name_lst = {}

    for name in names:
        mobj = re.match('\*\s*([a-zA-Q\.]+\d*\s+[a-zA-Z]{3}\s?[A-D]?)$', name)
        if mobj and name_parser.is_Bayer_name(mobj.group(1)):
            newname = name_parser.parse_Bayer_name(mobj.group(1))
            if 'Bayer' not in name_lst:
                name_lst['Bayer'] = []
            name_lst['Bayer'].append(newname)
            continue

        mobj = re.match('\*\s*(\d+\s+[a-zA-Z]{3}\s?[A-D]?)$', name)
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

        if name_parser.is_CPD_name(name):
            name_lst['CPD'] = name_parser.parse_CPD_name(name)
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

def select_main_name(name_lst, altname=None):
    if 'Bayer' in name_lst:
        main_name = name_lst['Bayer']
    elif 'Flamsteed' in name_lst:
        main_name = name_lst['Flamsteed']
    elif 'Var' in name_lst:
        main_name = name_lst['Var']
    elif 'HD' in name_lst:
        main_name = 'HD {}'.format(name_lst['HD'])
    elif 'DM' in name_lst:
        main_name = name_lst['DM']
    elif 'HIP' in name_lst:
        main_name = 'HIP {}'.format(name_lst['HIP'])
    elif 'CPD' in name_lst:
        main_name = name_lst['CPD']
    elif 'TYC' in name_lst:
        main_name = 'TYC {}'.format(name_lst['TYC'])
    elif 'GSC' in name_lst:
        main_name = 'GSC {}'.format(name_lst['GSC'])
    elif 'TIC' in name_lst:
        main_name = 'TIC {}'.format(name_lst['TIC'])
    elif '2MASS' in name_lst:
        main_name = '2MASS J{}'.format(name_lst['2MASS'])
    elif altname is not None:
        main_name = altname
    else:
        main_name = ''

    if isinstance(main_name, list):
        return main_name[0]
    else:
        return main_name

def make_crosstable(name_lst, coord_lst,
        columns=['HD','DM','HIP','TYC','TIC','GSC','CCDM','UCAC4',
                '2MASS','Gaia2'],
        verbose = False,
        ):

    dtype = [
        ('name',    'S50'),
        ('ra',      'f8'),
        ('dec',     'f8'),
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

    for i, (name, coord) in enumerate(zip(name_lst, coord_lst)):

        # query object ID to get names in various catalogues
        result1 = Simbad.query_objectids(name)
        if result1 is None:
            # cannot find this object in Simbad
            print('\033[31mWarning: Missing object id:{}\033[0m'.format(name))
            main_name = name
            mask_name = False
            if coord is None:
                ra,  mask_ra  = 0.0, True
                dec, mask_dec = 0.0, True
            else:
                ra,  mask_ra  = coord.ra.deg,  False
                dec, mask_dec = coord.dec.deg, False
            crossname_lst = get_names([main_name])
        else:
            crossname_lst = get_names(list(result1['ID']))

            # query object to get coordinates
            result2 = Simbad.query_object(name)
            if result2 is None:
                print('\033[31mWarning: Missing object:{}\033[0m'.format(name))
            row0 = result2[0]
            newcoord = SkyCoord(row0['RA'], row0['DEC'], unit=(units.hourangle, units.deg))

            # query additional catalogues to get more names
            # query TIC
            if 'TIC' in crossname_lst:
                tic = int(crossname_lst['TIC'])
                if verbose:
                    print(i, 'search TIC:', tic)
                result3 = Vizier(catalog='IV/38/tic', columns=['**'],
                                column_filters={'TIC':'={}'.format(tic)}
                                ).query_constraints()
                if len(result3)==0:
                    print('Error: no reults in TIC')
                    ticrow = None
                else:
                    ticrow = result3[0][0]
            else:
                if verbose:
                    print(i, 'search coordinate in TIC:',
                            newcoord.ra.deg, newcoord.dec.deg)
                result3 = Vizier(catalog='IV/38/tic', columns=['**','+_r']
                            ).query_region(newcoord, radius=3*units.arcsec)
                if len(result3)==0:
                    print('Error: no reults in TIC')
                    ticrow = None
                elif 'TYC' in crossname_lst:
                    mask = result3[0]['TYC']==crossname_lst['TYC']
                    if mask.sum()==1:
                        ticrow = result3[0][mask][0]
                    else:
                        ticrow = None
                elif len(result3[0])==1 and result3[0][0]['_r']<0.1:
                    print(result3[0])
                    ticrow = result3[0][0]
                else:
                    ticrow = None
            if ticrow is not None:
                if 'TIC' not in crossname_lst:
                    crossname_lst['TIC'] = ticrow['TIC']
                if 'HIP' not in crossname_lst and ticrow['HIP'] is not np.ma.masked:
                    crossname_lst['HIP'] = ticrow['HIP']
                if 'TYC' not in crossname_lst and ticrow['TYC'] is not np.ma.masked \
                    and len(ticrow['TYC'])>0:
                    crossname_lst['TYC'] = ticrow['TYC']
                if 'UCAC4' not in crossname_lst and ticrow['UCAC4'] is not np.ma.masked:
                    crossname_lst['UCAC4'] = ticrow['UCAC4']
                if '2MASS' not in crossname_lst and ticrow['_2MASS'] is not np.ma.masked:
                    crossname_lst['2MASS'] = ticrow['_2MASS']
                if 'Gaia2' not in crossname_lst and ticrow['GAIA'] is not np.ma.masked:
                    crossname_lst['Gaia2'] = ticrow['GAIA']
            
            main_name = select_main_name(crossname_lst, altname=name)
            mask_name = False
            ra,  mask_ra  = newcoord.ra.deg,  False
            dec, mask_dec = newcoord.dec.deg, False

        item = [
            (main_name, mask_name),
            (ra,        mask_ra),
            (dec,       mask_dec),
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

    crosstable['ra'].info.format  = '%9.5f'
    crosstable['dec'].info.format = '%9.5f'
    for column in ['DM','CCDM']:
        if column in crosstable.colnames:
            lenlist = [len(row[column]) for row in crosstable
                        if row[column] is not np.ma.masked]
            maxlen = max(lenlist)
            crosstable[column].info.format = '<{}s'.format(maxlen)

    return crosstable
