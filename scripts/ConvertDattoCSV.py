#!/usr/bin/env python

import numpy as np
import argparse
import os
import glob
import csv

# Script converts dat to csv files used by Python version of HAverager

########################
# setup input parameters
########################
parser = argparse.ArgumentParser(description="Converts dat to csv files used by Python version of HAverager")
parser.add_argument('fpath', metavar='fpath', type=str, nargs='+',
                    help='input file path.')

parser.add_argument("--Percent", action="store_true", dest="percent", default=False,
                    help="Presentation of input uncertainty. " +
                         "True: relative uncertainty in percent False: Absolute uncertainty")

arguments = parser.parse_args()


fpath = arguments.fpath
inpercent = arguments.percent

if len(fpath)==1:
    fnames = glob.glob(fpath[0])
else:
    fnames = fpath

# Loop over files
# Write output to .dat file
for fname in fnames:

    # Read the file
    print(fname)
    F = open(fname)

    nData = 0
    nCol = 0
    ColType = []
    ColName = []
    persent = ''
    counter = 0

    for line in F:
        counter+=1
        if ('NData' in line):
            nData = int(line.split('=')[1])

        if ('NColumn' in line):
            nCol = int(line.split('=')[1])

        if ('ColumnType' in line):
            l = (line.split('=')[1]).replace("'","").replace("\n","").split(',')
            for s in l:
                nl = s.split('*')
                if(len(nl)==2):
                    ColType.extend([nl[1]]*int(nl[0]))
                else:
                    ColType.append(s)
        if ('ColumnName' in line):
            ColInfo = (line.split('=')[1]).replace("'","").replace("\n","").split(',')
        if ('Percent' in line):
            s = (line.split('=')[1]).replace("'","").split(',')
            persent = [bool(x) for x in s]


        if ('END' in line):
            break

    data = np.loadtxt(fname, skiprows=counter)

    for i,s in enumerate(ColType):
        if 'Bin' in s:
            ColName.append('bin'+ColInfo[i].strip())
            continue
        if 'Sigma' in s:
            ColName.append('data')
            continue
        if 'Dummy' in s:
            data = np.delete(data, (i), axis=1)
            continue
        if 'uncor' in ColInfo[i]:
            ColName.append('stat')
        else:
            ColName.append(ColInfo[i].strip())


    # Write csv file
    f = open(os.path.basename(fname).replace('.dat', '.csv'), 'w')
    writer = csv.writer(f, delimiter=',', quoting=csv.QUOTE_NONE)
    writer.writerow(ColName)
    writer.writerows(data.astype('string'))
    f.close()
