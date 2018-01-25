import os
from ..utils.fitsio import get_bintable_info
from ..utils.asciitable import structitem_to_dict
from .base import _get_KIC_number

def find_KIC(name, output='dict'):
    '''
    Find records in Kepler Input Catalog (`V/133
    <http://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=V/133>`_, Kepler
    Mission Team, 2009).

    The data file used in this function is complied from the 10th version of
    KIC. It contains 13,161,029 records with consecutive KIC numbers. Proper
    motions are available for 12,944,973 objects, or 98% of the entire sample.
    Parallaxes are provided for 958 objects, and physical parameters
    (*T*:sub:`eff`, log\ *g*, log\ *Z* and *R*) are available for 2,106,821
    objects, or 16% of the entire sample.

    For more details, see :ref:`Kepler Input Catalog<catalog_kic>`.

    .. csv-table:: Descriptions of returned parameters
        :header: Key, Type, Unit, Description
        :widths: 30, 30, 30, 120

        KIC,    integer32, ,       KIC number
        RAdeg,  float64,   deg,    Right ascension (*α*) at J2000
        DEdeg,  float64,   deg,    Declination (*δ*) at J2000
        pmRA,   float32,   mas/yr, Proper motion in Right ascension with cos(*δ*) factor
        pmDE,   float32,   mas/yr, Proper motion in Declination
        Plx,    float32,   mas,    Parallax
        umag,   float32,   mag,    *u* magnitude in SDSS system
        gmag,   float32,   mag,    *g* magnitude in SDSS system
        rmag,   float32,   mag,    *r* magnitude in SDSS system
        imag,   float32,   mag,    *i* magnitude in SDSS system
        zmag,   float32,   mag,    *z* magnitude in SDSS system
        grmag,  float32,   mag,    Magnitude in GRed band
        d51mag, float32,   mag,    Magnitude in DDO-51 filter
        kepmag, float32,   mag,    Magnitude in Kepler band
        flag_g, integer16, ,       "Galaxy flag (0 for star, 1 for galaxy)"
        flag_v, integer16, ,       "Variable flag (0 for normal, 1 for variable)"
        cq,     string5,   ,       Origin of Kepelr magnitude
        fv,     integer16, ,       "0 for outside Kepler FOV, 1/2 for inside, 2 for Kepler target"
        Teff,   integer16, K,      Effective temperature
        logg,   float32,   dex,    Surface gravity
        FeH,    float32,   dex,    Metallicity
        EBV,    float32,   mag,    Color excess in *B* − *V*
        Av,     float32,   mag,    Extinction in *V* magnitude
        R,      float32,   Rsun,   Stellar radius

    Args:
        name (string or integer): Name or number of star.
        output (string): Type of output results. Either *"dict"* or *"dtype"*
            (:class:`numpy.dtype`).
    Returns:
        dict or :class:`numpy.dtype`: Record in catalogue.
    Raises:
        FileNotExist: Catalogue file does not exist
        UnrecognizedName: Input name can not be recognized
    Examples:
        Find *K*:sub:`p` magnitude of Kepler-13 (KOI-13, KIC 9941662)

        .. code-block:: python
        
            from stella.catalog.find_catalog import find_KIC
    
            res = find_KIC('KIC 9941662')
            print(res['kepmag'])
            # output:
            # 9.958000183105469

    '''

    filename = os.path.join(os.getenv('STELLA_DATA'), 'catalog/KIC.fits')
    if not os.path.exists(filename):
        raise FileNotExist(filename)

    nbyte, nrow, ncol, pos, dtype, fmtfunc = get_bintable_info(filename)
    infile = open(filename,'rb')

    kic = _get_KIC_number(name)
    if kic>0 and kic<=nrow:
        infile.seek(pos+(kic-1)*nbyte,0)
        item = fmtfunc(infile.read(nbyte))
    else:
        pass

    infile.close()

    if output == 'ndarray':
        return item
    elif output == 'dict':
        return structitem_to_dict(item)
    else:
        return None

