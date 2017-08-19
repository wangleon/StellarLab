from __future__ import print_function
import os
import math
import struct
import numpy as np
import astropy.io.fits as fits

from ..utils.fitsio import get_bintable_info
from ..utils.asciitable import structitem_to_dict
from .errors import FileNotExist, ItemNotFound, UnrecognizedName
from .name import _get_HIP_number, _get_KIC_number


def _search_HIP_catalogue(names, filename, epoch, output):
    '''
    Search data in either HIP catalogue (`I/239
    <http://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=I/239>`_, Perryman+
    1997) or HIP New Reduction (`I/311
    <http://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=I/311>`_, van Leeuwen
    2007).

    Parameters
    -----------
    names : *string*, *integer*, *list*, *tuple*, or *numpy.ndarray*
        Names of stars
    filename : *string*
        Name of catalogue file. Either *"HIP.fits"* or *"HIP2.fits"*
    epoch : *float*
        Epoch of astrometric parameters. Default is 2000.0
    output : *'dict'* (default) or *'ndarray'*
        Data type of output. Default is *'dict'*

    Returns
    --------
    row : *numpy.ndarray* or *dict*
        Record array or dictionary
    '''

    # check table file
    catfile = os.path.join(os.getenv('STELLA_DATA'), 'catalog/%s'%filename)
    if not os.path.exists(catfile):
        raise FileNotExist(catfile)

    nbyte, nrow, ncol, pos, dtype, fmtfunc = get_bintable_info(catfile)
    infile = open(catfile,'rb')

    def change_epoch(item, epoch):
        pm_ra = item['pmRA']*1e-3/3600. # convert pm_RA from mas/yr to deg/yr
        pm_de = item['pmDE']*1e-3/3600. # convert pm_DE from mas/yr to deg/yr
        item['RAdeg'] += (epoch-1991.25)*pm_ra/math.cos(item['DEdeg']/180.*math.pi)
        item['DEdeg'] += (epoch-1991.25)*pm_de

    if isinstance(names, list) or \
       isinstance(names, tuple) or \
       (isinstance(names, np.ndarray) and np.issubdtype(names.dtype, int)):

        hip_lst = list(map(_get_HIP_number, names))

        result_lst = []
        for hip in hip_lst:
            if hip is None:
                # hip = 672 is the common null record in both HIP and HIP New
                infile.seek(pos+(672-1)*nbyte,0)
                item = fmtfunc(infile.read(nbyte))
            else:
                infile.seek(pos+(hip-1)*nbyte,0)
                item = fmtfunc(infile.read(nbyte))
                change_epoch(item, epoch)
            result_lst.append(item)

        infile.close()

        if output == 'ndarray':
            return np.array(result_lst, dtype=item.dtype)
        elif output == 'dict':
            return list(map(structitem_to_dict, result_lst))
        else:
            return None

    else:
        # input is a single starname

        hip = _get_HIP_number(names)
        if hip is None:
            # hip = 672 is the common null record in both HIP and HIP New
            infile.seek(pos+(672-1)*nbyte,0)
            item = fmtfunc(infile.read(nbyte))
        else:
            infile.seek(pos+(hip-1)*nbyte,0)
            item = fmtfunc(infile.read(nbyte))
            change_epoch(item, epoch)

        infile.close()

        if output == 'ndarray':
            return item
        elif output == 'dict':
            return structitem_to_dict(item)
        else:
            return None

def find_HIP(names, epoch=2000.0, output='dict'):
    '''
    Find records in Hipparcos Catalogue (`I/239
    <http://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=I/239>`_, Perryman+
    1997).

    Parameters
    -----------
    names : *string*, *integer*, *list*, *tuple*, or *numpy.ndarray*
        Names of stars
    epoch : *float*
        Epoch of astrometric parameters. Default is 2000.0
    output : *'dict'* (default) or *'ndarray'*
        Data type of output. Default is *'dict'*

    Returns
    --------
    row : *numpy.ndarray* or *dict*
        Record array or dictionary
    
    Raises
    -------
    FileNotExist
        Catalogue file does not exist

    Notes
    ------
    The Hipparcos Catalogue (I/239, Perryman et al. 1997) is the product of the
    Hipparcos satellite of the European Space Agency (ESA). During its four
    years of operation from Novermber 1989 to March 1993, the satellite measured
    accurate triangular parallaxes and proper motions for a large number of
    stars. The published astrometric solutions are given in International
    Celestial Reference System (ICRS) at epoch J1991.25.

    This catalogue contains 118,218 records, with HIP numbers ranging from 1 to
    120,416. The limiting magnitude is *V* ~ 12.4, and the *V* magnitude
    completeness is 7.3 ~ 9.0, depending on the galactic latitude and spectral
    type. The median error in parallax is 0.97 mas for stars with Hipparcos
    magnitudes brighter than 9.

    .. csv-table:: Descriptions of returned parameters
        :header: "Key", "Type", "Unit", "Description"
        :widths: 30, 30, 30, 100

        "HIP",      "integer32", "",       "HIP number"
        "RAdeg",    "float64",   "deg",    "Right ascension (*α*) in ICRS at given epoch"
        "DEdeg",    "float64",   "deg",    "Declination (*δ*) in ICRS at at given epoch"
        "Vmag",     "float32",   "mag",    "*V* magnitude in Johnson system"
        "Plx",      "float32",   "mas",    "Parallax"
        "e_Plx",    "float32",   "mas",    "Error in parallax"
        "pmRA",     "float32",   "mas/yr", "Proper motion in Right ascension"
        "pmDE",     "float32",   "mas/yr", "Proper motion in Declination"
        "e_pmRA",   "float32",   "mas/yr", "Error in proper motion in RA"
        "e_pmDE",   "float32",   "mas/yr", "Error in proper motion in Dec"
        "BTmag",    "float32",   "mag",    "Mean *B*:sub:`T` magnitude"
        "e_BTmag",  "float32",   "mag",    "Error in *B*:sub:`T` magnitude"
        "VTmag",    "float32",   "mag",    "Mean *V*:sub:`T` magnitude"
        "e_VTmag",  "float32",   "mag",    "Error *V*:sub:`T` magnitude"
        "B-V",      "float32",   "mag",    "*B* − *V* color in Johnson system"
        "e_B-V",    "float32",   "mag",    "Error in *B* − *V* color"
        "r_B-V",    "character", "",       "Source of *B* − *V* color"
        "V-I",      "float32",   "mag",    "*V* − *I* color in Cousin system"
        "e_V-I",    "float32",   "mag",    "Error in *V* − *I* color"
        "r_V-I",    "character", "",       "Source of *V* − *I* color"
        "Hpmag",    "float32",   "mag",    "Median magnitude in Hipparcos system"
        "e_Hpmag",  "float32",   "mag",    "Error in *Hp* magnitude"
        "Hpscat",   "float32",   "mag",    "Scatter on *Hp* magnitude"
        "o_Hpmag",  "integer16", "",       "Number of observations for *Hp* magnitude"
        "SpType",   "string",    "mag",    "Spectral type"
        "r_SpType", "character", "",       "Source of Spectral type"

    Examples
    --------

    .. code-block:: python

        import numpy as np
        from stella.catalog.find_catalog import find_HIP

        # find the parametres for tau Cet (HIP 8102)
        item = find_HIP(8102, epoch=1991.25)

        # find the parameters for HIP 65000 ~ HIP 65399
        hip_lst = np.arange(65000, 65400)
        res = find_HIP(hip_lst, output='ndarray')


    References
    ----------
    * `Perryman et al., 1997, A&A, 323, L49 <http://adsabs.harvard.edu/abs/1997A&A...323L..49P>`_
    '''

    return _search_HIP_catalogue(names, 'HIP.fits', epoch, output)

def find_HIP2(names, epoch=2000.0, output='dict'):
    '''
    Find record in Hipparcos New Reduction (`I/311
    <http://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=I/311>`_, van Leeuwen
    2007).

    Parameters
    -----------
    names : *string*, *integer*, *list*, *tuple*, or *numpy.ndarray*
        Names of stars
    epoch : *float*
        Epoch of astrometric parameters. Default is 2000.0
    output : *'dict'* (default) or *'ndarray'*
        Data type of output. Default is *'dict'*

    Returns
    --------
    row : *numpy.ndarray* or *dict*
        Record array or dictionary

    Raises
    -------
    FileNotExist
        Catalogue file does not exist

    Notes
    ------
    This catalogue is the results of the new reduction (`van Leeuwen, 2007
    <http://adsabs.harvard.edu/abs/2007A&A...474..653V>`_) of the astrometric
    data obtained by the Hipparcos satellite.

    This catalogue contains 117,955 records, with HIP numbers ranging from 1 to
    120,404. The astrometric solution are provided in ICRS at epoch J1991.25.
    However, parameters at any epoch can be returned with the `epoch` argument.
    Default is *epoch* = 2000.0.
    
    .. csv-table:: Descriptions of returned parameters
        :header: "Key", "Type", "Unit", "Description"
        :widths: 30, 30, 30, 100

        "HIP",      "integer32", "",       "HIP number"
        "RAdeg",    "float64",   "deg",    "Right ascension (*α*) in ICRS at given epoch"
        "DEdeg",    "float64",   "deg",    "Declination (*δ*) in ICRS at at given epoch"
        "Plx",      "float32",   "mas",    "Parallax"
        "e_Plx",    "float32",   "mas",    "Error in parallax"
        "pmRA",     "float32",   "mas/yr", "Proper motion in Right ascension"
        "pmDE",     "float32",   "mas/yr", "Proper motion in Declination"
        "e_pmRA",   "float32",   "mas/yr", "Error in proper motion in RA"
        "e_pmDE",   "float32",   "mas/yr", "Error in proper motion in Dec"
        "B-V",      "float32",   "mag",    "*B* − *V* color in Johnson system"
        "e_B-V",    "float32",   "mag",    "Error in *B* − *V* color"
        "V-I",      "float32",   "mag",    "*V* − *I* color in Cousin system"
        "Hpmag",    "float32",   "mag",    "Median magnitude in Hipparcos system"
        "e_Hpmag",  "float32",   "mag",    "Error in *Hp* magnitude"
        "Hpscat",   "float32",   "mag",    "Scatter on *Hp* magnitude"

    Examples
    ---------
    Find the ICRS coordinate of τ Ceti (HIP 8102)

    .. code-block:: python
    
        from stella.catalog.find_catalog import find_HIP2

        res1 = find_HIP2(8102, epoch=1991.25)
        res2 = find_HIP2(8102)
        print(res1['RAdeg'], res1['DEdeg'])
        print(res2['RAdeg'], res2['DEdeg'])
        # output:
        # 26.021364586713265 -15.939555724635493
        # 26.017014215650022 -15.937479641367434



    References
    ----------
    * `van Leeuwen, 2007, A&A, 474, 653 <http://adsabs.harvard.edu/abs/2007A&A...474..653V>`_

    '''

    return _search_HIP_catalogue(names, 'HIP2.fits', epoch, output)

def find_TYC(starname):

    filename = os.path.join(os.getenv('STELLA_DATA'), 'catalog/TYC.fits')
    if not os.path.exists(filename):
        raise FileNotExist(filename)

    # get tyc1, tyc2, and tyc3
    if starname[0:3]=='TYC':
        g = starname[3:].split('-')
        tyc1 = int(g[0])
        tyc2 = int(g[1])
        tyc3 = int(g[2])

    # get information of FITS table
    nbyte, nrow, ncol, pos, dtype, fmtfunc = get_bintable_info(filename)

    target = np.int32((tyc1<<18) + (tyc2<<4) + (tyc3<<1))

    infile = open(filename)
    i1, i2 = 0, nrow-1
    while(i2-i1 > 1):
        i3 = int((i1+i2)/2)
        infile.seek(pos + i3*nbyte, 0)
        key = struct.unpack('>i',infile.read(4))[0]
        if target < key:
            i2 = i3
        elif target > key:
            i1 = i3
        else:
            break
    infile.seek(-4,1)
    row = fmtfunc(infile.read(nbyte))
    
    # looking for possible companion
    if i3 != nrow - 1:
        key2 = struct.unpack('>i',infile.read(4))[0]
        if key2 == key + 1:
            infile.seek(-4,1)
            row2 = fmtfunc(infile.read(nbyte))
            print('Warning: There are more than 1 star matched')
            t1 = (key2 & 0b11111111111111000000000000000000)/2**18
            t2 = (key2 & 0b00000000000000111111111111110000)/2**4
            t3 = (key2 & 0b00000000000000000000000000001110)/2
            print(t1,t2,t3,row2)
    infile.close()
    return row

def find_KIC(names, output='dict'):
    '''
    Find records in Kepler Input Catalog (`V/133
    <http://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=V/133>`_, Kepler
    Mission Team, 2009).

    Parameters
    -----------
    names : *string*, *integer*, *list*, *tuple*, or *numpy.ndarray*
        Names of stars
    output : *'dict'* (default) or *'ndarray'*
        Data type of output. Default is *'dict'*
        
    Returns
    --------
    row : *numpy.ndarray* or *dict*
        Record array or dictionary

    Raises
    -------
    FileNotExist
        Catalogue file does not exist
    UnrecognizedName
        Input name can not be recognized

    Notes
    -----
    The Kepler Input Catalog (V/133, Kepler Mission Team, 2009) contains
    photometric and physical data of over 13 million objects in the Kepler field
    of view. For details, see `Brown et al. 2011
    <http://adsabs.harvard.edu/abs/2011AJ....142..112B>`_.

    This catalogue is complied from the 10th version of KIC. It contains
    13,161,029 records with consecutive KIC numbers. Proper motions are
    available for 12,944,973 objects, or 98% of the entire sample. Parallaxes
    are provided for 958 objects, and physical parameters (*T*:sub:`eff`,
    log\ *g*, [Fe/H] and *R*) are available for 2,106,821 objects, or 16% of the
    entire sample.

    There are some known biases on the physical parameters of the Kepler Input
    Catalog. We refer the readers to `Pinsonneault et al. 2012
    <http://adsabs.harvard.edu/abs/2012ApJS..199...30P>`_ for an improved
    temperature scale and `Dong et al. 2014
    <http://adsabs.harvard.edu/abs/2014ApJ...789L...3D>`_ for an improved
    metallicity scale. Besides, `Huber et al. 2014
    <http://adsabs.harvard.edu/abs/2014ApJS..211....2H>`_ derived parameters for
    196,468 stars using Kepler photometric data, including 11,532 unclassified
    targes in this catalogue.


    .. csv-table:: Descriptions of returned parameters
        :header: "Key", "Type", "Unit", "Description"
        :widths: 30, 30, 30, 120

        "KIC",    "integer32", "",       "KIC number"
        "RAdeg",  "float64",   "deg",    "Right ascension (*α*) at J2000"
        "DEdeg",  "float64",   "deg",    "Declination (*δ*) at J2000"
        "pmRA",   "float32",   "mas/yr", "Proper motion in Right ascension"
        "pmDE",   "float32",   "mas/yr", "Proper motion in Declination"
        "Plx",    "float32",   "mas",    "Parallax"
        "umag",   "float32",   "mag",    "*u* magnitude in SDSS system"
        "gmag",   "float32",   "mag",    "*g* magnitude in SDSS system"
        "rmag",   "float32",   "mag",    "*r* magnitude in SDSS system"
        "imag",   "float32",   "mag",    "*i* magnitude in SDSS system"
        "zmag",   "float32",   "mag",    "*z* magnitude in SDSS system"
        "grmag",  "float32",   "mag",    "Magnitude in GRed band"
        "d51mag", "float32",   "mag",    "Magnitude in DDO-51 filter"
        "kepmag", "float32",   "mag",    "Magnitude in Kepler band"
        "flag_g", "integer16", "",       "Galaxy flag (0 for star, 1 for galaxy)"
        "flag_v", "integer16", "",       "Variable flag (0 for normal, 1 for variable)"
        "cq",     "string5",   "",       "Origin of Kepelr magnitude"
        "fv",     "integer16", "",       "0 for outside Kepler FOV, 1/2 for inside, 2 for Kepler target"
        "Teff",   "integer16", "K",      "Effective temperature"
        "logg",   "float32",   "dex",    "Surface gravity"
        "FeH",    "float32",   "dex",    "Metallicity"
        "EBV",    "float32",   "mag",    "Color excess in *B* − *V*"
        "Av",     "float32",   "mag",    "Extinction in *V* magnitude"
        "R",      "float32",   "Rsun",   "Stellar radius"

    References
    -----------
    * `Brown et al., 2011, AJ, 142, 112 <http://adsabs.harvard.edu/abs/2011AJ....142..112B>`_
    * `Dong et al., 2014, ApJ, 789, L3 <http://adsabs.harvard.edu/abs/2014ApJ...789L...3D>`_
    * `Huber et al., 2014, ApJS, 211, 2 <http://adsabs.harvard.edu/abs/2014ApJS..211....2H>`_
    * `Pinsonneault et al., 2012, ApJS, 199, 30 <http://adsabs.harvard.edu/abs/2012ApJS..199...30P>`_
    '''

    filename = os.path.join(os.getenv('STELLA_DATA'), 'catalog/KIC.fits')
    if not os.path.exists(filename):
        raise FileNotExist(filename)

    nbyte, nrow, ncol, pos, dtype, fmtfunc = get_bintable_info(filename)
    infile = open(filename,'rb')

    if isinstance(names, list) or \
       isinstance(names, tuple) or \
       (isinstance(names, np.ndarray) and np.issubdtype(names.dtype, int)):

        kic_lst = list(map(_get_KIC_number, names))

        result_lst = []
        for kic in kic_lst:
            if kic>0 and kic<=nrow:
                infile.seek(pos+(kic-1)*nbyte,0)
                item = fmtfunc(infile.read(nbyte))
                result_lst.append(item)
            else:
                pass

        if output == 'ndarray':
            return np.array(result_lst, dtype=item.dtype)
        elif output == 'dict':
            return list(map(structitem_to_dict, result_lst))
        else:
            return None

    else:
        # input is a single starname

        kic = _get_KIC_number(names)
        if kic>0 and kic<=nrow:
            infile.seek(pos+(kic-1)*nbyte,0)
            item = fmtfunc(infile.read(nbyte))
        else:
            pass

        if output == 'ndarray':
            return item
        elif output == 'dict':
            return structitem_to_dict(item)
        else:
            return None
