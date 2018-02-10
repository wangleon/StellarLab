import os
import math
import numpy as np
import astropy.io.fits as fits
from ..utils.fitsio import get_bintable_info
from ..utils.asciitable import structitem_to_dict
from .base import _get_HIP_number

catfile1 = os.path.join(os.getenv('STELLA_DATA'), 'catalog/HIP.fits')
catfile2 = os.path.join(os.getenv('STELLA_DATA'), 'catalog/HIP2.fits')

def _find_HIP_object(name, catfile, data_info, epoch=2000.0, output='dict'):
    '''
    Find record for an object in either HIP catalogue (`I/239
    <http://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=I/239>`_, Perryman+
    1997) or HIP New Reduction (`I/311
    <http://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=I/311>`_, van Leeuwen
    2007).

    Args:
        name (string or integer): Name or number of star.
        catfile (string): Name of the catalogue file.
        data_info (dict): Information of the FITS catalogue file.
        epoch (float): Epoch of the output astrometric parameters.
        output (string): Type of output results. Either *"dict"* or *"dtype"*
            (:class:`numpy.dtype`).
    Returns:
        dict or :class:`numpy.dtype`: Record in catalogue.
    '''

    def change_epoch(item, epoch):
        pm_ra = item['pmRA']*1e-3/3600. # convert pm_RA from mas/yr to deg/yr
        pm_de = item['pmDE']*1e-3/3600. # convert pm_DE from mas/yr to deg/yr
        item['RAdeg'] += (epoch-1991.25)*pm_ra/math.cos(item['DEdeg']/180.*math.pi)
        item['DEdeg'] += (epoch-1991.25)*pm_de

    hip = _get_HIP_number(name)

    pos     = data_info['pos']
    nbyte   = data_info['nbyte']
    fmtfunc = data_info['fmtfunc']

    infile = open(catfile,'rb')

    if hip is None:
        # return a null result
        # hip = 672 is the common null record in both HIP and HIP New
        infile.seek(pos+(672-1)*nbyte, 0)
        item = fmtfunc(infile.read(nbyte))
    else:
        infile.seek(pos+(hip-1)*nbyte, 0)
        item = fmtfunc(infile.read(nbyte))
        change_epoch(item, epoch)
    infile.close()

    if output == 'dtype':
        return item
    elif output == 'dict':
        return structitem_to_dict(item)
    else:
        return None


class _HIP(object):
    '''Class for Hipparcos Catalogue (`I/239
    <http://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=I/239>`_, Perryman+
    1997).

    The HIP catalogue contains 118,218 records, with HIP numbers ranging
    from 1 to 120,416.
    For more details, see :ref:`Hipparcos Catalogue<catalog_hip>`.
    '''

    def __init__(self):
        self.catfile = catfile1
        self._data_info = None

    def _get_data_info(self):
        '''Get information of FITS table.'''
        nbyte, nrow, ncol, pos, dtype, fmtfunc = get_bintable_info(self.catfile)
        self._data_info = {
                'nbyte'  : nbyte,
                'nrow'   : nrow,
                'ncol'   : ncol,
                'pos'    : pos,
                'dtype'  : dtype,
                'fmtfunc': fmtfunc
                }

    def find_object(self, name, epoch=2000.0, output='dict'):
        '''
        Find record for an object in Hipparcos Catalogue.
    
        .. csv-table:: Descriptions of returned parameters
            :header: Key, Type, Unit, Description
            :widths: 30, 30, 30, 100
    
            HIP,      integer32, ,       HIP number
            RAdeg,    float64,   deg,    Right ascension (*α*) in ICRS at given epoch
            DEdeg,    float64,   deg,    Declination (*δ*) in ICRS at at given epoch
            Vmag,     float32,   mag,    *V* magnitude in Johnson system
            Plx,      float32,   mas,    Parallax
            e_Plx,    float32,   mas,    Error in parallax
            pmRA,     float32,   mas/yr, Proper motion in Right ascension with cos(*δ*) factor
            pmDE,     float32,   mas/yr, Proper motion in Declination
            e_pmRA,   float32,   mas/yr, Error in proper motion in RA
            e_pmDE,   float32,   mas/yr, Error in proper motion in Dec
            BTmag,    float32,   mag,    Mean *B*:sub:`T` magnitude
            e_BTmag,  float32,   mag,    Error in *B*:sub:`T` magnitude
            VTmag,    float32,   mag,    Mean *V*:sub:`T` magnitude
            e_VTmag,  float32,   mag,    Error in *V*:sub:`T` magnitude
            B-V,      float32,   mag,    *B* − *V* color in Johnson system
            e_B-V,    float32,   mag,    Error in *B* − *V* color
            r_B-V,    character, ,       Source of *B* − *V* color
            V-I,      float32,   mag,    *V* − *I* color in Cousin system
            e_V-I,    float32,   mag,    Error in *V* − *I* color
            r_V-I,    character, ,       Source of *V* − *I* color
            Hpmag,    float32,   mag,    Median magnitude in Hipparcos system
            e_Hpmag,  float32,   mag,    Error in *Hp* magnitude
            Hpscat,   float32,   mag,    Scatter on *Hp* magnitude
            o_Hpmag,  integer16, ,       Number of observations for *Hp* magnitude
            SpType,   string,    mag,    Spectral type
            r_SpType, character, ,       Source of Spectral type
    
        Args:
            name (string or integer): Name or number of star.
            epoch (float): Epoch of output astrometric parameters.
            output (string): Type of output results. Either *"dict"* or
                *"dtype"* (:class:`numpy.dtype`).
        Returns:
            dict or :class:`numpy.dtype`: Record in catalogue.
        Raises: 
            FileNotExist: Catalogue file does not exist
        Examples:
            Find record for τ Ceti (HIP 8102)
    
            .. code-block:: python
    
                >>> from stella.catalog import HIP
                >>> rec = HIP.find_object(8102, epoch=1991.25)
                >>> rec['Vmag'], rec['B-V'], rec['Plx'], rec['e_Plx']
                (3.490000009536743, 0.7269999980926514, 274.1700134277344, 0.800000011920929)
    
        '''

        if self._data_info is None:
            self._get_data_info()

        return _find_HIP_object(name, self.catfile, self._data_info, epoch, output)


class _HIP2(object):
    '''Class for Hipparcos New Reduction (`I/311
    <http://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=I/311>`_, van Leeuwen
    2007).

    The Hipparcos Catalogue New Reduction contains 117,955 records, with HIP
    numbers ranging from 1 to 120,404.
    For more details, see :ref:`Hipparcos Catalogue New Reduction<catalog_hip2>`.
    '''

    def __init__(self):
        self.catfile = catfile2
        self._data_info = None

    def _get_data_info(self):
        '''Get information of FITS table.'''
        nbyte, nrow, ncol, pos, dtype, fmtfunc = get_bintable_info(self.catfile)
        self._data_info = {
                'nbyte'  : nbyte,
                'nrow'   : nrow,
                'ncol'   : ncol,
                'pos'    : pos,
                'dtype'  : dtype,
                'fmtfunc': fmtfunc
                }

    def find_object(self, name, epoch=2000.0, output='dict'):
        '''
        Find record for an object in Hipparcos Catalogue New Reduction.

        .. csv-table:: Descriptions of returned parameters
            :header: Key, Type, Unit, Description
            :widths: 30, 30, 30, 100
    
            HIP,     integer32, ,       HIP number
            RAdeg,   float64,   deg,    Right ascension (*α*) in ICRS at given epoch
            DEdeg,   float64,   deg,    Declination (*δ*) in ICRS at at given epoch
            Plx,     float32,   mas,    Parallax
            e_Plx,   float32,   mas,    Error in parallax
            pmRA,    float32,   mas/yr, Proper motion in Right ascension with cos(*δ*) factor
            pmDE,    float32,   mas/yr, Proper motion in Declination
            e_pmRA,  float32,   mas/yr, Error in proper motion in RA
            e_pmDE,  float32,   mas/yr, Error in proper motion in Dec
            B-V,     float32,   mag,    *B* − *V* color in Johnson system
            e_B-V,   float32,   mag,    Error in *B* − *V* color
            V-I,     float32,   mag,    *V* − *I* color in Cousin system
            Hpmag,   float32,   mag,    Median magnitude in Hipparcos system
            e_Hpmag, float32,   mag,    Error in *Hp* magnitude
            Hpscat,  float32,   mag,    Scatter on *Hp* magnitude
    
        Args:
            name (string or integer): Name or number of star.
            epoch (float): Epoch of output astrometric parameters.
            output (string): Type of output results. Either *"dict"* or *"dtype"*
                (:class:`numpy.dtype`).
        Returns:
            dict or :class:`numpy.dtype`: Record in catalogue.
        Raises:
            FileNotExist: Catalogue file does not exist
        Examples:
            Find the ICRS coordinate of τ Ceti (HIP 8102)
    
            .. code-block:: python
            
                >>> from stella.catalog import HIP2
                >>> rec1 = HIP2.find_object(8102, epoch=1991.25)
                >>> rec2 = HIP2.find_object(8102)
                >>> rec1['RAdeg'], rec1['DEdeg']
                (26.021364586713265, -15.939555724635493)
                >>> rec2['RAdeg'], rec2['DEdeg']
                (26.017014215650022, -15.937479641367434)

        '''

        if self._data_info is None:
            self._get_data_info()

        return _find_HIP_object(name, self.catfile, self._data_info, epoch, output)

HIP = _HIP()
HIP2 = _HIP2()

def convert_HIP_to_fits(inputfile, outputfile='HIP.fits'):
    '''Convert ASCII catalog files of HIP catalogue to FITS table.

    Args:
        inputfile (string): Name of the input ASCII file
        outputfile (string): Name of the output FITS file
    Returns:
        No returns.

    Examples:

        .. code-block:: python

            >>> import os
            >>> from stella.catalog.hip import convert_HIP_to_fits
            >>> inputfile = os.path.join(os.getenv('ASTRO_DATA'), 'catalog/I/239/hip_main.dat')
            >>> convert_HIP_to_fits(inputfile, 'HIP.fits')
    '''
    from .base import _str_to_float

    types = [
            ('HIP',    np.int32),
            ('RAdeg',  np.float64),
            ('DEdeg',  np.float64),
            ('Vmag',   np.float32),
            ('Plx',    np.float32),
            ('e_Plx',  np.float32),
            ('pmRA',   np.float32),
            ('pmDE',   np.float32),
            ('e_pmRA', np.float32),
            ('e_pmDE', np.float32),
            ('BTmag',  np.float32),
            ('e_BTmag',np.float32),
            ('VTmag',  np.float32),
            ('e_VTmag',np.float32),
            ('B-V',    np.float32),
            ('e_B-V',  np.float32),
            ('r_B-V',  'S1'),
            ('V-I',    np.float32),
            ('e_V-I',  np.float32),
            ('r_V-I',  'S1'),
            ('Hpmag',  np.float32),
            ('e_Hpmag',np.float32),
            ('Hpscat', np.float32),
            ('o_Hpmag',np.int16),
            #('CCDM',   'S10'),
            #('HD',     np.int32),
            #('BD',     'S10'),
            #('CoD',    'S10'),
            ('SpType', 'S12'),
            ('r_SpType','S1'),
            ]
    tmp = list(zip(*types))
    record = np.dtype({'names':tmp[0],'formats':tmp[1]})

    fill_item = np.array((0, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN,
                             np.NaN, np.NaN, np.NaN, np.NaN, np.NaN,
                             np.NaN, np.NaN, np.NaN, np.NaN, np.NaN,
                             '',     np.NaN, np.NaN, '', np.NaN, np.NaN, np.NaN,
                             -32768, '',''),dtype=record)

    data = {}

    infile = open(inputfile)
    for row in infile:
        hip = int(row[8:14])
        rah = int(row[17:19])
        ram = int(row[20:22])
        ras = float(row[23:28])
        radeg = (rah + ram/60. + ras/3600.)*15.
        ded = abs(int(row[30:32]))
        dem = int(row[33:35])
        des = float(row[36:40])
        dedeg = ded + dem/60. + des/3600.
        if row[29]=='-':
            dedeg = -dedeg
        vmag     = _str_to_float(row[41:46], np.NaN)
        plx      = _str_to_float(row[79:86], np.NaN)
        e_plx    = _str_to_float(row[119:125], np.NaN)
        pmRA     = _str_to_float(row[87:95], np.NaN)
        pmDE     = _str_to_float(row[96:104], np.NaN)
        e_pmRA   = _str_to_float(row[126:132], np.NaN)
        e_pmDE   = _str_to_float(row[133:139], np.NaN)
        BTmag    = _str_to_float(row[217:223], np.NaN)
        VTmag    = _str_to_float(row[230:236], np.NaN)
        e_BTmag  = _str_to_float(row[224:229], np.NaN)
        e_VTmag  = _str_to_float(row[237:242], np.NaN)
        BV       = _str_to_float(row[245:251], np.NaN)
        e_BV     = _str_to_float(row[252:257], np.NaN)
        r_BV     = row[258].strip()
        VI       = _str_to_float(row[260:264], np.NaN)
        e_VI     = _str_to_float(row[265:269], np.NaN)
        r_VI     = row[270].strip()
        Hpmag    = _str_to_float(row[274:281], np.NaN)
        e_Hpmag  = _str_to_float(row[282:288], np.NaN)
        Hpscat   = _str_to_float(row[289:294], np.NaN)
        if row[295:298].strip()=='':
            o_Hpmag = 0
        else:
            o_Hpmag = int(row[295:298])

        if not np.isnan(pmRA):
            pm_ra = pmRA*1e-3/3600. # convert pm_RA  from mas/yr to degree/yr
            #radeg += (2000.0-1991.25)*pm_ra/math.cos(dedeg/180.*math.pi)

        if not np.isnan(pmDE):
            pm_de = pmDE*1e-3/3600. # convert pm_Dec from mas/yr to degree/yr
            #dedeg += (2000.0-1991.25)*pm_de

        SpType   = row[435:447].strip()
        r_SpType = row[448].strip()

        item = np.array((hip, radeg, dedeg, vmag, plx, e_plx,
                         pmRA, pmDE, e_pmRA, e_pmDE,
                         BTmag, e_BTmag, VTmag, e_VTmag,
                         BV, e_BV, r_BV, VI, e_VI, r_VI,
                         Hpmag, e_Hpmag, Hpscat, o_Hpmag,
                         #CCDM, HD, BD, CoD,
                         SpType, r_SpType,
                        ),dtype=record)
        if hip in data:
            print('Error: Duplicate Records for HIP', hip)
        data[hip] = item
    infile.close()

    newdata = []
    for hip in range(1, max(data.keys())+1):
        if hip in data:
            newdata.append(data[hip])
        else:
            newdata.append(fill_item)
    newdata = np.array(newdata, dtype=record)

    pri_hdu = fits.PrimaryHDU()
    tbl_hdu = fits.BinTableHDU(newdata)
    hdu_lst = fits.HDUList([pri_hdu,tbl_hdu])

    if os.path.exists(outputfile):
        os.remove(outputfile)
    hdu_lst.writeto(outputfile)

def convert_HIP2_to_fits(inputfile, outputfile='HIP2.fits'):
    '''Convert ASCII catalog files of HIP catalogue (New Reductionto) FITS
    table.

    Args:
        inputfile (string): Name of the input ASCII file
        outputfile (string): Name of the output FITS file
    Returns:
        No returns.

    Examples:

        .. code-block:: python

            >>> import os
            >>> from stella.catalog.hip import convert_HIP2_to_fits
            >>> inputfile = os.path.join(os.getenv('ASTRO_DATA'), 'catalog/I/311/hip2.dat')
            >>> convert_HIP2_to_fits(inputfile, 'HIP2.fits')
    '''
    from .base import _str_to_float
    types = [
            ('HIP',    np.int32),
            ('RAdeg',  np.float64),
            ('DEdeg',  np.float64),
            ('Plx',    np.float32),
            ('e_Plx',  np.float32),
            ('pmRA',   np.float32),
            ('pmDE',   np.float32),
            ('e_pmRA', np.float32),
            ('e_pmDE', np.float32),
            ('B-V',    np.float32),
            ('e_B-V',  np.float32),
            ('V-I',    np.float32),
            ('Hpmag',  np.float32),
            ('e_Hpmag',np.float32),
            ('Hpscat', np.float32),
            ]
    tmp = list(zip(*types))
    record = np.dtype({'names':tmp[0],'formats':tmp[1]})

    fill_item = np.array(
                    (0, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN,
                        np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN),
                    dtype=record)

    data = {}

    infile = open(inputfile)
    for row in infile:
        hip     = int(row[0:6])
        radeg   = float(row[15:28])/math.pi*180.
        dedeg   = float(row[29:42])/math.pi*180.
        plx     = _str_to_float(row[43:50], np.NaN)
        e_plx   = _str_to_float(row[83:89], np.NaN)
        pmRA    = _str_to_float(row[51:59], np.NaN)
        pmDE    = _str_to_float(row[60:68], np.NaN)
        e_pmRA  = _str_to_float(row[90:96], np.NaN)
        e_pmDE  = _str_to_float(row[97:103], np.NaN)
        Hpmag   = _str_to_float(row[129:136], np.NaN)
        e_Hpmag = _str_to_float(row[137:143], np.NaN)
        Hpscat  = _str_to_float(row[144:149], np.NaN)
        BV      = _str_to_float(row[152:158], np.NaN)
        e_BV    = _str_to_float(row[159:164], np.NaN)
        VI      = _str_to_float(row[165:171], np.NaN)

        if not np.isnan(pmRA):
            pm_ra = pmRA*1e-3/3600. # convert pm_RA  from mas/yr to degree/yr
            #radeg += (2000.0-1991.25)*pm_ra/math.cos(dedeg/180.*math.pi)

        if not np.isnan(pmDE):
            pm_de = pmDE*1e-3/3600. # convert pm_Dec from mas/yr to degree/yr
            #dedeg += (2000.0-1991.25)*pm_de

        item = np.array((hip,radeg, dedeg, plx, e_plx,
                        pmRA, pmDE, e_pmRA, e_pmDE,
                        BV, e_BV, VI, Hpmag,e_Hpmag, Hpscat), dtype=record)
        if hip in data:
            print('Error: Duplicate record for HIP', hip)
        data[hip] = item

    infile.close()

    count_good, count_null = 0, 0
    newdata = []
    for hip in range(1, max(data.keys())+1):
        if hip in data:
            newdata.append(data[hip])
            count_good += 1
        else:
            print(hip)
            newdata.append(fill_item)
            count_null += 1
    newdata = np.array(newdata, dtype=record)

    pri_hdu = fits.PrimaryHDU()
    tbl_hdu = fits.BinTableHDU(newdata)
    hdu_lst = fits.HDUList([pri_hdu,tbl_hdu])

    if os.path.exists(outputfile):
        os.remove(outputfile)
    hdu_lst.writeto(outputfile)

