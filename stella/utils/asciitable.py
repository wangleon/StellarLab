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
    add_names, add_formats, add_values = [], [], []

    for row in file1:
        row = row.strip()
        if len(row) == 0 or row[0] == '#':
            continue

        if row[0] == '%':
            row = row[1:].strip()
            if '=' in row:
                g = row.split('=')
                key   = g[0].strip()
                value = g[1].strip()
                if key == 'names':
                    names = [v.strip() for v in value.split(',')]
                elif key == 'formats':
                    formats = [v.strip() for v in value.split(',')]
                elif key == 'delimiter':
                    delimiter = value.strip()
                elif key[0:6] == 'global':
                    g2 = key.split()
                    name = g2[2]
                    fmt = g2[1]
                    if fmt == 'string':
                        value = value
                        fmt = 'S%d'%len(value)
                    elif fmt == 'float':
                        value = float(value)
                        fmt = 'f'
                    elif fmt == 'double':
                        value = float(value)
                        fmt = 'd'
                    elif fmt == 'int':
                        if len(value)>2 and value[0:2]=='0b':
                            value = int(value, 2)
                        elif len(value)>2 and value[0:2]=='0x':
                            value = int(value, 16)
                        elif len(value)>1 and value[0]=='0':
                            value = int(value, 8)
                        else:
                            value = int(value)
                        fmt = 'i'
                    else:
                        raise ValueError
                    add_names.append(name)
                    add_formats.append(fmt)
                    add_values.append(value)
                else:
                    pass
        else:
            # data begin
            # initialize dytpe
            if custom is None and None not in [names, formats]:
                ncol = len(names)
                if len(add_names)>0:
                    for name, fmt in zip(add_names, add_formats):
                        names.append(name)
                        formats.append(fmt)
                custom = np.dtype({'names':list(names),'formats':list(formats)})
                data = []

            if delimiter is None:
                g = row.split()
            else:
                g = row.split(delimiter)
            rec = []
            for i, descr in enumerate(custom.descr[0:ncol]):
                value = g[i].strip()
                if descr[1][1]=='i':
                    if len(value)>2 and value[0:2]=='0b':
                        value = int(value, 2)
                    elif len(value)>2 and value[0:2]=='0x':
                        value = int(value, 16)
                    elif len(value)>1 and value[0]=='0':
                        value = int(value, 8)
                    else:
                        value = int(value)
                elif descr[1][1] in ['f','d']:
                    try:
                        value = float(value)
                    except:
                        value = None
                elif descr[1][1]=='S':
                    value = value
                else:
                    print 'Unknow describer: ',descr[1]
                rec.append(value)
            for value in add_values:
                rec.append(value)
            recdata = np.array(tuple(rec), dtype=custom)
            data.append(recdata)
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

