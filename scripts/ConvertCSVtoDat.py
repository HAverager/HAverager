#!/usr/bin/env python

import numpy as np
import pandas as pd
from optparse import OptionParser
import os
import glob

########################
# setup input parameters
########################
usage = ("DatasetGen.py parameters \nCreate dataset with gaussian uncertainties.  \n")
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

    print(fname)
    df = pd.read_csv(fname)

    bins = df.columns.values
    vals = df.values

    nData = vals.shape[0]
    nSyst = vals.shape[1] - 3

    f = open(os.path.basename(fname).replace('.csv', '.dat'), 'w')

    f.write('&Data\n')
    f.write('   Name = \'Data%i\'\n' % nData)
    f.write('   NData = %i\n' % nData)
    f.write('   NColumn = %i\n' % (nSyst + 3))
    f.write('   ColumnType = ')
    for b in bins:
        if 'bin' in b:
            f.write('\'Bin\', ')
        elif 'data' in b:
            f.write('\'Sigma\', ')
        else:
            f.write('\'Error\', ')
    f.write('\n')
    f.write('   ColumnName = ')
    for s in bins:
        f.write("'{}', ".format(s))
    f.write('\n')
    f.write('   Reaction = \'Bla\'\n')
    f.write('   Percent = true')
    for s in range(nSyst):
        f.write(',true')
    f.write('\n')
    f.write('&END\n')

    # Loop over data point
    for d in range(nData):
        f.write('%4.0f ' % vals[d, 0])
        f.write('%15.10f ' % vals[d, 1])
        f.write('%8.3f ' % vals[d, 2])
        # Loop over systematics
        for s in range(nSyst):
            f.write('%8.3f ' % vals[d, s + 3])
        f.write('\n')
    f.write('\n')
    f.close()
