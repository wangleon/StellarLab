import re

constellations = {
        'And': ['Andromedae',           'Andromeda'             ],
        'Ant': ['Antliae',              'Antlia'                ],
        'Aps': ['Apodis',               'Apus'                  ],
        'Aqr': ['Aquarii',              'Aquarius'              ],
        'Aql': ['Aquilae',              'Aquila'                ],
        'Ara': ['Arae',                 'Ara'                   ],
        'Ari': ['Arietis',              'Aries'                 ],
        'Aur': ['Aurigae',              'Auriga'                ],
        'Boo': ['Bootis',               'Bootes'                ],
        'Cae': ['Caeli',                'Caelum'                ],
        'Cam': ['Camelopardalis',       'Camelopardalis'        ],
        'Cnc': ['Cancri',               'Cancer'                ],
        'CVn': ['Canum Venaticorum',    'Canes Venatici'        ],
        'CMa': ['Canis Majoris',        'Canis Major'           ],
        'CMi': ['Canis Minoris',        'Canis Minor'           ],
        'Cap': ['Capricorni',           'Capricornus'           ],
        'Car': ['Carinae',              'Carina'                ],
        'Cas': ['Cassiopeiae',          'Cassiopeia'            ],
        'Cen': ['Centauri',             'Centaurus'             ],
        'Cep': ['Cephei',               'Cepheus'               ],
        'Cet': ['Ceti',                 'Cetus'                 ],
        'Cha': ['Chamaeleotis',         'Chamaeleon'            ],
        'Cir': ['Circini',              'Circinus'              ],
        'Col': ['Columbae',             'Columba'               ],
        'Com': ['Comae Berenices',      'Coma Berenices'        ],
        'CrA': ['Coronae Australis',    'Corona Australis'      ],
        'CrB': ['Coronae Borealis',     'Corona Borealis'       ],
        'Crv': ['Corvi',                'Corvus'                ],
        'Crt': ['Crateris',             'Crater'                ],
        'Cru': ['Crucis',               'Crux'                  ],
        'Cyg': ['Cygni',                'Cygnus'                ],
        'Del': ['Delphini',             'Delphinus'             ],
        'Dor': ['Doradus',              'Dorado'                ],
        'Dra': ['Draconis',             'Draco'                 ],
        'Equ': ['Equulei',              'Equuleus'              ],
        'Eri': ['Eridani',              'Eridanus'              ],
        'For': ['Fornacis',             'Fornax'                ],
        'Gem': ['Geminorum',            'Gemini'                ],
        'Gru': ['Gruis',                'Grus'                  ],
        'Her': ['Herculis',             'Hercules'              ],
        'Hor': ['Horologii',            'Horologium'            ],
        'Hya': ['Hydrae',               'Hydra'                 ],
        'Hyi': ['Hydri',                'Hydrus'                ],
        'Ind': ['Indi',                 'Indus'                 ],
        'Lac': ['Lacertae',             'Lacerta'               ],
        'Leo': ['Leonis',               'Leo'                   ],
        'LMi': ['Leonis Minoris',       'Leo Minor'             ],
        'Lep': ['Leporus',              'Lepus'                 ],
        'Lib': ['Librae',               'Libra'                 ],
        'Lup': ['Lupi',                 'Lupus'                 ],
        'Lyn': ['Lyncis',               'Lynx'                  ],
        'Lyr': ['Lyrae',                'Lyra'                  ],
        'Men': ['Mensae',               'Mensa'                 ],
        'Mic': ['Microscopii',          'Microscopium'          ],
        'Mon': ['Monocerotis',          'Monoceros'             ],
        'Mus': ['Muscae',               'Musca'                 ],
        'Nor': ['Normae',               'Norma'                 ],
        'Oct': ['Octantis',             'Octans'                ],
        'Oph': ['Ophiuchi',             'Ophiuchus'             ],
        'Ori': ['Orionis',              'Orion'                 ],
        'Pav': ['Pavonis',              'Pavo'                  ],
        'Peg': ['Pegasi',               'Pegasus'               ],
        'Per': ['Persei',               'Perseus'               ],
        'Phe': ['Phoenicis',            'Phoenix'               ],
        'Pic': ['Pictoris',             'Pictor'                ],
        'Psc': ['Piscium',              'Pisces'                ],
        'PsA': ['Piscis Austrini',      'Piscis Austrinus'      ],
        'Pup': ['Puppis',               'Puppis'                ],
        'Pyx': ['Pyxidis',              'Pyxis'                 ],
        'Ret': ['Reticuli',             'Reticulum'             ],
        'Sge': ['Sagittae',             'Sagitta'               ],
        'Sgr': ['Sagittarii',           'Sagittarius'           ],
        'Sco': ['Scorpii',              'Scorpius'              ],
        'Scl': ['Sculptoris',           'Sculptor'              ],
        'Sct': ['Scuti',                'Scutum'                ],
        'Ser': ['Serpentis',            'Serpens'               ],
        'Sex': ['Sextantis',            'Sextans'               ],
        'Tau': ['Tauri',                'Taurus'                ],
        'Tel': ['Telescopii',           'Telescopium'           ],
        'Tri': ['Trianguli',            'Triangulum'            ],
        'TrA': ['Trianguli Australis',  'Triangulum Australe'   ],
        'Tuc': ['Tucanae',              'Tucana'                ],
        'UMa': ['Ursae Majoris',        'Ursa Major'            ],
        'UMi': ['Ursae Minoris',        'Ursa Minor'            ],
        'Vel': ['Velorum',              'Vela'                  ],
        'Vir': ['Virginis',             'Virgo'                 ],
        'Vol': ['Volantis',             'Volans'                ],
        'Vul': ['Vulpeculae',           'Vulpecula'             ],
        }

greek_letters = {
        'alf': 'alpha',
        'bet': 'beta',
        'gam': 'gamma',
        'det': 'delta',
        'eps': 'epsilon',
        'zet': 'zeta',
        'eta': 'eta',
        'the': 'theta',
        'iot': 'iota',
        'kap': 'kappa',
        'lam': 'lambda',
        'mu' : 'mu',
        'nu' : 'nu',
        'xi' : 'xi',
        'omi': 'omicron',
        'pi' : 'pi',
        'rho': 'rho',
        'sig': 'sigma',
        'tau': 'tau',
        'ups': 'upsilon',
        'phi': 'phi',
        'chi': 'chi',
        'psi': 'psi',
        'ome': 'omega',
        }

def _get_star_number1(name, key):
    '''
    Convert star name with the form of `SSS NNNN` to its integer number `NNNN`.

    Args:
        name (integer or string): Name of a star.
        key (string): Prefix of the star name.
    Returns:
        integer: Number of the star in the catalog. If fail a *None* value will
            be returned.

    '''
    if isinstance(name, int):
        return name
    elif isinstance(name, str):
        if name[0:len(key)] == key:
            return int(name[len(key):])
        elif name.isdigit():
            return int(name)
        else:
            return None
    else:
        return None

def _get_HIP_number(name):
    '''Convert star name in *Hipparcos Catalogue* to an integer HIP number.

    Args:
        name (string or integer): Name of the star.
    Returns:
        integer: HIP number.
    '''
    return _get_star_number1(name, 'HIP')

def _get_KIC_number(name):
    '''Convert star name in *Kepler Input Catalog* to an integer KIC number.

    Args:
        name (string or integer): Name of the star
    Returns:
        integer: KIC number
    '''
    return _get_star_number1(name, 'KIC')

def _get_EPIC_number(name):
    '''Convert star name in *K2 Ecliptic Plane Input Catalog* to an integer EPIC
    number.

    Args:
        name (string or integer): Name of the star
    Returns:
        integer: EPIC number
    '''
    return _get_star_number1(name, 'EPIC')

def _get_TYC_number(name):
    '''Convert star name in *Tycho-2 Catalogue* to TYC numbers (TYC1, TYC2,
    TYC3).

    Args:
        name (string): Name of the star.
    Returns:
        tuple: A tuple of TYC numbers (TYC1, TYC2, TYC3).
    '''
    if name[0:3]=='TYC':
        g = name[3:].split('-')
        tyc1, tyc2, tyc3 = int(g[0]), int(g[1]), int(g[2])
        return (tyc1, tyc2, tyc3)
    elif len(name.split('-'))==3:
        g = name.split('-')
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

def _get_regular_HIP_name(name):
    '''Convert an HIP name in *Hipparcos Catalogue* to its regular form
    `"HIP NNN"`.

    Args:
        name (string or integer): Name or HIP number of a star (e.g. `"HIP8276"`,
            `"HIP 8276"`, `8443`).
    Returns:
        string: Regular HD name `"HIP NNNN"`.
    See also:
        * :ref:`catalog_hip`
    '''

    if isinstance(name, str):
        name = name.strip()
        return 'HIP ' + name[3:].strip()
    elif isinstance(name, int):
        return 'HIP %d'%name
    else:
        raise ValueError

def _get_regular_HD_name(name):
    '''Convert an HD name in *Henry Draper Catalogue* to its regular form
    `"HD NNNN"` or `"HD NNNN C"`.
    
    Args:
        name (string or integer): Name or HD number of a star (e.g. `"HD8276"`,
            `"HD 8276A"`, `8443`).
    Returns:
        string: Regular HD name `"HD NNNN"` or `"HD NNNN C"`.
    See also:

        * `Henry Draper Catalogue and Extension (III/125A)
          <http://vizier.u-strasbg.fr/cgi-bin/VizieR?-source=III/135A/catalog>`_

    '''
    if isinstance(name, str):
        name = name.strip()
        if name[-1].isalpha():
            comp = name[-1]
            return 'HD %d %s'%(int(name[2:-1]), comp)
        else:
            return 'HD %d'%(int(name[2:]))
    elif isinstance(name, int):
        return 'HD %d'%name
    else:
        raise ValueError


def _get_regular_BD_name(name):
    '''Convert a BD name in *Bonner Durchmusterung* to its regular form 
    `"BD+MM NNNN"`.
    '''
    name = name.strip()
    # parse 'BD+33 23' to 'BD +33 23'
    name = 'BD '+ name[2:].strip()
    # parse 'BD 33 23' to 'BD +33 23'
    g = name.split()
    if g[1][0].isdigit():
        g[1] = '+'+g[1]
        name = ' '.join(g)

    # parse 'BD +33 23A' to 'BD +33 23 A'
    # but keep 'BD +33 23a'
    if name[-1].isalpha():
        comp = name[-1]
        if comp.isupper():
            name = name[:-1].strip()+' '+comp
        elif comp.islower():
            name = name[:-1].strip()+comp
        else:
            raise ValueError

    # parse 'BD -3 23' to 'BD -03 23'
    g = name.split()
    if len(g[1])!=3:
        pm = g[1][0]
        num = abs(int(g[1]))
        g[1] = pm + str(num).rjust(2,'0')

    # parse 'BD -03 0023' to 'BD -03 23'
    g[2] = str(int(g[2]))
    name = ' '.join(g)

    return name

def _get_regular_CD_name(name):
    '''Convert a CD name in *Cordoba Durchmusterung* to its regular form
    `"CD+MM NNNN"`.
    '''

    name = name.strip()

    # parse 'CD-33 23' to 'CD -33 23'
    name = 'CD '+name[2:].strip()

    # parse 'CD -33 23A' to 'CD -33 23 A'
    if name[-1].isalpha():
        comp = name[-1]
        if comp.isupper():
            name = name[:-1].strip()+' '+comp
        elif comp.islower():
            name = name[:-1].strip()+comp
        else:
            raise ValueError

    # parse 'CD -22 0023' to 'CD -22 23'
    g = name.split()
    g[2] = str(int(g[2]))
    name = ' '.join(g)

    return name


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
    '''Convert a TYC name in *Tycho-2 Catalogue* to its regular form
    `"TYC NNN-NNN-N"`.
    '''
    if len(args)==1 and isinstance(args[0], str):
        name = args[0].strip()
        return 'TYC '+name[3:].strip()
    elif len(args)==3:
        return 'TYC '+'-'.join([str(args[0]),str(args[1]),str(args[2])])
    else:
        return None
