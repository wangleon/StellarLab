#!/usr/bin/env python3
import os

path = os.path.join(os.getenv('ASTRO_DATA'), 'catalog/V/133/kic/')

filename_lst = []
for direct in 'sn':
    for i in range(90):
        fn = os.path.join(path, '%s%02d.dat'%(direct,i))
        if os.path.exists(fn):
            filename_lst.append(fn)
if len(filename_lst)==0:
    print('Error: Cannot find catalog file in %s'%path)

count_pm = 0
count_plx = 0
count_para = 0
for filename in filename_lst:
    kic_lst = []
    infile = open(filename)
    for row in infile:
        if row[0]== '#':
            continue
        kic = int(row[0:8])
        kic_lst.append(kic)
        if len(row[30:38].strip())!=0 or len(row[39:47].strip())!=0:
            count_pm += 1
        if len(row[48:56].strip())!=0:
            count_plx += 1
        if len(row)>225:
            count_para += 1
    infile.close()
    #print(os.path.basename(filename), min(kic_lst), max(kic_lst))
print('N(PM)=%d, N(Plx)=%d, N(para)=%d'%(count_pm, count_plx, count_para))
