#!/usr/bin/env python

import numpy as np
import pandas as pd
from optparse import OptionParser
import os
import glob
import csv

# Script converts csv to dat files used by fotran version of HAverager

########################
# setup input parameters
########################
usage = ("ConvertCSVtoDat.py parameters \nConvert csv to dat files used by fotran version of HAverager \n")
parser = OptionParser(usage)
parser.add_option("-f", dest="fpath", type="string", help="input file path")
parser.add_option("--InPercent", "-i", action="store_true", dest="inpercent", default=False,
                  help="Presentation of input uncertainty. "
                       + "True: relative uncertainty in %, False: Absolute uncertainty")

options, arguments = parser.parse_args()

fpath = options.fpath
inpercent = options.inpercent

fnames = glob.glob(fpath)

# Loop over files
# Write output to .dat file
for fname in fnames:

    # Read the file
    print fname
    F = open(fname)

    nData = 0
    nCol = 0
    ColType = []
    ColName = ''
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
            ColName = (line.split('=')[1]).replace("'","").replace("\n","").split(',')
        if ('Percent' in line):
            s = (line.split('=')[1]).replace("'","").split(',')
            persent = [bool(x) for x in s]


        if ('END' in line):
            break


    for i,s in enumerate(ColType):
        if 'Bin' in s:
            ColName[i] = 'bin'+ColName[i].strip()
        if 'Sigma' in s:
            ColName[i] = 'data'
        if 'uncor' in ColName[i]:
            ColName[i] = 'stat'

    data = np.loadtxt(fname, skiprows=counter)
    print data

    # Write csv file
    f = open(os.path.basename(fname).replace('.dat', '.csv'), 'w')
    writer = csv.writer(f, delimiter=',', quoting=csv.QUOTE_NONE)
    writer.writerow(ColName)
    writer.writerows(data.astype('string'))
    f.close()