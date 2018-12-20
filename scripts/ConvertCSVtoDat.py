#!/usr/bin/env python

import pandas as pd
import argparse
import os
import glob

# Script converts csv to dat files used by fotran version of HAverager

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
        f.write('%8.3f ' % vals[d, 0])
        f.write('%15.10f ' % vals[d, 1])
        f.write('%8.3f ' % vals[d, 2])
        # Loop over systematics
        for s in range(nSyst):
            f.write('%8.3f ' % vals[d, s + 3])
        f.write('\n')
    f.write('\n')
    f.close()
