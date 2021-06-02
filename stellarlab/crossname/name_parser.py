import re

constellation_abbr_lst = [
    'And', 'Ant', 'Aps', 'Aqr', 'Aql', 'Ara', 'Ari', 'Aur', 'Boo', 'Cae',
    'Cam', 'Cnc', 'CVn', 'CMa', 'CMi', 'Cap', 'Car', 'Cas', 'Cen', 'Cep',
    'Cet', 'Cha', 'Cir', 'Col', 'Com', 'CrA', 'CrB', 'Crv', 'Crt', 'Cru',
    'Cyg', 'Del', 'Dor', 'Dra', 'Equ', 'Eri', 'For', 'Gem', 'Gru', 'Her',
    'Hor', 'Hya', 'Hyi', 'Ind', 'Lac', 'Leo', 'LMi', 'Lep', 'Lib', 'Lup',
    'Lyn', 'Lyr', 'Men', 'Mic', 'Mon', 'Mus', 'Nor', 'Oct', 'Oph', 'Ori',
    'Pav', 'Peg', 'Per', 'Phe', 'Pic', 'Psc', 'PsA', 'Pup', 'Pyx', 'Ret',
    'Sge', 'Sgr', 'Sco', 'Scl', 'Sct', 'Ser', 'Sex', 'Tau', 'Tel', 'Tri',
    'TrA', 'Tuc', 'UMa', 'UMi', 'Vel', 'Vir', 'Vol', 'Vul',
    ]


greek_letter_abbr_lst = {
    'alf': 'alpha',     'bet': 'beta',      'gam': 'gamma',
    'del': 'delta',     'eps': 'epsilon',   'zet': 'zeta',
    'eta': 'eta',       'the': 'theta',     'tet': 'theta',
    'iot': 'iota',      'kap': 'kappa',     'lam': 'lambda',
    'mu.': 'mu',        'nu.': 'nu',        'ksi': 'xi',
    'omi': 'omicron',   'pi.': 'pi',        'rho': 'rho',
    'sig': 'sigma',     'tau': 'tau',       'ups': 'upsilon',
    'phi': 'phi',       'chi': 'chi',       'psi': 'psi',
    'ome': 'omega',
    }

greek_letter_lst = greek_letter_abbr_lst.values()

def is_Bayer_name(starname):
    mobj = re.match('^([a-zA-Q\.]+)(\d*)\s*([a-zA-Z]{3})$', starname)
    if mobj:
        letter = mobj.group(1)
        number = mobj.group(2)
        conste = mobj.group(3)
        if conste in constellation_abbr_lst and \
            (letter in greek_letter_abbr_lst \
            or letter in greek_letter_lst \
            or re.match('^[a-zA-Q]$', letter)):
            return True
        else:
            return False
    return False

def parse_Bayer_name(starname):
    mobj = re.match('^([a-zA-Q\.]+)(\d*)\s*([a-zA-Z]{3})$', starname)
    if mobj:
        letter = mobj.group(1)
        number = mobj.group(2)
        conste = mobj.group(3)
        if letter in greek_letter_abbr_lst:
            letter = greek_letter_abbr_lst[letter]
        return '{}{} {}'.format(letter, number, conste)
    else:
        return ''

def is_Flamsteed_name(starname):

    mobj = re.match('^(\d+)\s+([a-zA-Z]{3})$', starname)
    if mobj:
        number = mobj.group(1)
        conste = mobj.group(2)
        if conste in constellation_abbr_lst:
            return True
        else:
            return False
    return False

def parse_Flamsteed_name(starname):
    
    mobj = re.match('^(\d+)\s+([a-zA-Z]{3})$', starname)
    if mobj:
        number = mobj.group(1)
        conste = mobj.group(2)
        return '{} {}'.format(number, conste)
    else:
        return ''

def is_DM_name(starname):
    if re.match('([BC]D[+\-]?\d+\s*\d+[a-zA-Z]?)$', starname):
        return True
    else:
        return False

def parse_DM_name(starname):
    mobj = re.match('([BC]D)([+\-]?\d+)\s*(\d+)([a-zA-Z]?)$', starname)
    if mobj:
        catalog = mobj.group(1)
        zone    = int(mobj.group(2))
        number  = int(mobj.group(3))
        comp    = mobj.group(4)
        return '{}{:=+03d} {:d}{}'.format(catalog, zone, number, comp)
    else:
        return ''

def is_Var_name(starname):
    mobj = re.match('([R-Z])\s*([a-zA-Z]{3})$', starname)
    if mobj:
        # variable names from Q to Z
        number = mobj.group(1)
        conste = mobj.group(2)
        if conste in constellation_abbr_lst:
            return True
        else:
            return False

    mobj = re.match('([A-Z]{2})\s*([a-zA-Z]{3})$', starname)
    if mobj:
        # Variable names with double uppercase letters
        number = mobj.group(1)
        conste = mobj.group(2)
        if conste in constellation_abbr_lst:
            return True
        else:
            return False

    mobj = re.match('(V\d+)\s*([a-zA-Z]{3})$', starname)
    if mobj:
        # Variable names V???
        number = mobj.group(1)
        conste = mobj.group(2)
        if conste in constellation_abbr_lst:
            return True
        else:
            return False

    return False

def parse_Var_name(starname):
    mobj = re.match('([R-Z])\s*([a-zA-Z]{3})$', starname)
    if mobj:
        # variable names from Q to Z
        number = mobj.group(1)
        conste = mobj.group(2)
        if conste in constellation_abbr_lst:
            return '{} {}'.format(number, conste)
        else:
            return ''

    mobj = re.match('([A-Z]{2})\s*([a-zA-Z]{3})$', starname)
    if mobj:
        # Variable names with double uppercase letters
        number = mobj.group(1)
        conste = mobj.group(2)
        if conste in constellation_abbr_lst:
            return '{} {}'.format(number, conste)
        else:
            return ''

    mobj =  re.match('(V\d+)\s*([a-zA-Z]{3})$', starname)
    if mobj:
        # Variable names V???
        number = mobj.group(1)
        conste = mobj.group(2)
        if conste in constellation_abbr_lst:
            return '{} {}'.format(number, conste)
        else:
            return ''

    return ''


def parse_TYC_name(starname):
    mobj = re.match('TYC\s*(\d*)\-(\d*)\-(\d)$', starname)
    if mobj:
        tyc1 = int(mobj.group(1))
        tyc2 = int(mobj.group(2))
        tyc3 = int(mobj.group(3))
        return '{:04d}-{:05d}-{:d}'.format(tyc1, tyc2, tyc3)
    else:
        return ''
