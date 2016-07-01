import struct
import numpy as np

def tform_to_format(tform):
    '''
    Convert TFORM string in FITS binary table to fmt in Python struct.
    '''
    if tform == 'L': return 'b' # 1 byte, boolean
    if tform == 'B': return 'B' # 1 byte, unsigned byte
    if tform == 'I': return 'h' # 2 bytes, integer
    if tform == 'J': return 'i' # 4 bytes, integer
    if tform == 'K': return 'l' # 8 bytes, integer
    if tform == 'E': return 'f' # 4 bytes, single-precision float
    if tform == 'D': return 'd' # 8 bytes, double-precision float
    if tform == 'C': return '?' # 8 bytes, single-precision complex
    if tform == 'M': return '?' # 16 bytes, double-precision complex
    if tform == 'P': return '?' # 8 bytes, 32 bits array descriptor
    if tform == 'Q': return '?' # 16 bytes, 64 bits array descriptor
    if tform[-1] == 'A': return tform[0:-1]+'s' # 1 bytes, character


def tform_to_dtype(tform):
    '''
    Convert TFORM string in FITS binary table to format Python array.
    '''
    if tform == 'L': return np.bool # 1 byte, boolean
    #if tform == 'B': return 'B' # 1 byte, unsigned byte
    if tform == 'I': return np.int16 # 2 bytes, integer
    if tform == 'J': return np.int32 # 4 bytes, integer
    if tform == 'K': return np.int64 # 8 bytes, integer
    if tform == 'E': return np.float32 # 4 bytes, single-precision float
    if tform == 'D': return np.float64 # 8 bytes, double-precision float
    if tform == 'C': return '?' # 8 bytes, single-precision complex
    if tform == 'M': return '?' # 16 bytes, double-precision complex
    if tform == 'P': return '?' # 8 bytes, 32 bits array descriptor
    if tform == 'Q': return '?' # 16 bytes, 64 bits array descriptor
    if tform[-1] == 'A': return 'S'+tform[0:-1] # 1 bytes, character


def get_bintable_info(filename,extension=1):
    '''
    Return the information of the binary table in FITS file.
    Parameters
    ----------
    filename : string
        Name of the FITS file
    extension : integer
        Which extension is the binary table in
    Returns
    -------
    naxis1 : integer
    naxis2 : integer
    tfields : list
    position : integer
    dtype : 
    fmtfunc : function
    '''
    infile = open(filename)
    current_hdu = 0
    while(True):
        block = infile.read(36*80)
        if block[0:8] == 'XTENSION':
            current_hdu += 1
            if block[10:30].strip()[1:-1]=='BINTABLE' and \
                current_hdu==extension:
                infile.seek(-36*80,1)
                break
    count = 0
    while(True):
        row = infile.read(80)
        count += 1
        if row[0:3]=='END':
            infile.seek((36-count%36)*80,1)
            break
        elif row[0:6]=='NAXIS1':
            naxis1 = int(row[10:30])
        elif row[0:6]=='NAXIS2':
            naxis2 = int(row[10:30])
        elif row[0:7]=='TFIELDS':
            tfields = int(row[10:30])
            ttype_lst = ['' for j in range(tfields)]
            tform_lst = ['' for j in range(tfields)]
        elif row[0:5]=='TTYPE':
            index = int(row[5:8])
            ttype_lst[index-1] = row[10:30].strip()[1:-1].strip()
        elif row[0:5]=='TFORM':
            index = int(row[5:8])
            tform_lst[index-1] = row[10:30].strip()[1:-1].strip()

    position = infile.tell()
    infile.close()

    formats = tuple([tform_to_dtype(v) for v in tform_lst])
    dtype = np.dtype({'names':ttype_lst,'formats':formats})

    fmt = '>'+(''.join([tform_to_format(v) for v in tform_lst]))
    fmtfunc = lambda string: np.array(struct.unpack(fmt, string),dtype=dtype)
    return naxis1,naxis2,tfields,position,dtype,fmtfunc

