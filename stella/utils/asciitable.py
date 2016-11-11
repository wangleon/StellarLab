#!/usr/bin/env python
import os
import numpy as np

def load_txt(filename):
    '''
    Load structured array from an ascii file.
    '''

    # open file
    file1 = open(filename)

    # initialize parameters
    names, formats, custom, delimiter, data = None, None, None, None, None

    for row in file1:
        row = row.strip()

        if len(row)==0 or row[0]=='#':
            # comment line or empty line
            continue

        elif row[0]=='%':

            if '=' not in row:
                continue
            key   = row[1:].split('=')[0].strip()
            value = row[1:].split('=')[1].strip()
            if key == 'names':
                names = [v.strip() for v in value.split(',')]
            elif key == 'formats':
                formats = [v.strip() for v in value.split(',')]
            elif key == 'delimiter':
                delimiter = value.strip()
            else:
                pass
            continue

        # if names and formats are ready, initialize dytpe
        # execute only once.
        if custom is None and names is not None and formats is not None:
            custom = np.dtype({'names':list(names), 'formats': list(formats)})
            data = []

        if custom is not None:
            if delimiter is None:
                g = row.split()
            else:
                g = row.split(delimiter)
            rec = []
            for i, descr in enumerate(custom.descr):
                if descr[1][1]=='i':
                    try:
                        d = int(g[i])
                    except:
                        d = None
                elif descr[1][1] in ['f','d']:
                    try:
                        d = float(g[i])
                    except:
                        d = None
                elif descr[1][1]=='S':
                    d = g[i]
                else:
                    print 'Unknow describer: ',descr[1]
                rec.append(d)
            recdata = np.array(tuple(rec), dtype=custom)
            data.append(recdata)
        else:
            print 'Numpy.dtype not defined'
    file1.close()
    data = np.array(data, dtype=custom)
    return data

def save_txt(filename, data, format_dict=None):
    '''
    Save numpy structured array to an ascii file
    '''

    file1 = open(filename, 'w')

    names, formats = [], []
    for rec in data.dtype.descr:
        names.append(rec[0])
        formats.append(rec[1])

    # find dlimiter
    for delimiter in [' ','|',',']:
        has_delimiter = False
        for rec in data.dtype.descr:
            if rec[1][1]=='S':
                for s in data[:][rec[0]]:
                    if delimiter in s:
                        has_delimiter = True
                        break
            if has_delimiter:
                break
        if has_delimiter:
            continue
        else:
            break
    if has_delimiter:
        print "can't find seperator"
        raise ValueError

    # write delimiter, names and formats.
    if delimiter != ' ':
        file1.write('%% delimiter = %s%s'%(delimiter, os.linesep))
    file1.write('%% names = %s%s'%(', '.join(names), os.linesep))
    file1.write('%% formats = %s%s'%(', '.join(formats), os.linesep))

    # find max lengths of each column
    maxlens = [0 for v in names]
    for row in data:
        for i,s in enumerate(row):
            if len(str(s))>maxlens[i]:
                maxlens[i] = len(str(s))

    # write data
    if format_dict == None:
        for row in data:
            lst = [str(v).ljust(maxlens[i]) for i,v in enumerate(list(row))]
            file1.write('%s%s'%(delimiter.join(lst), os.linesep))
    else:
        for row in data:
            lst = []
            for i,v in enumerate(list(row)):
                if names[i] in format_dict:
                    string = format_dict[names[i]]%v
                else:
                    string = str(v).ljust(maxlens[i])
                lst.append(string)
            file1.write('%s%s'%(delimiter.join(lst), os.linesep))

    file1.close()

