#!/usr/bin/env python

import os
import time
import matplotlib.pyplot as plt

# Import data reader and averager
import averager
import DataReader

# Script to test time performance of HAvergaer with different options

# perform averaging
def RunAverager(UseBlas=False, nToyMC=0, SysImp=False, itr=0):

    # read the data
    data, stat, syst = DataReader.paverage('', '', '')

    # initialization (optional information)
    averager.avin.cleaninvars()

    averager.avin.initvariables()
    averager.avin.setoutputfolder('./TOutP')
    averager.avin.initeration = itr
    averager.avin.inwriteoriginal = True

    averager.avin.infixstat = False
    averager.avin.incorrectstatbias = False
    averager.avin.inrescalestatsep = False

    averager.avin.inpostrotatesyst = False
    averager.avin.indosystimpact = SysImp
    averager.avin.inntoymc = nToyMC
    averager.avin.inuseblas = UseBlas

    start = time.time()
    averager.average(data, stat, syst)
    end = time.time()
    return (end-start)

# Test 1.
nsyst = [20, 50, 100, 150, 200, 250, 300, 450]
t1 = []
t2 = []

for ns in nsyst:
    os.system('DatasetGen.py -d 700 -s {}'.format(ns))
    t1.append(RunAverager(UseBlas=True))
    t2.append(RunAverager())


plt.figure()
plt.plot(nsyst, t1, label='with blas')
plt.plot(nsyst, t2, label='without blas')
plt.legend(numpoints=1, loc=0)
plt.xlabel('number of systematics')
plt.ylabel('running time')
plt.savefig('TimevsSyst.pdf')

# Test 2.
ndata = [100, 200, 500, 1000, 3000]
t1 = []
t2 = []

for ns in ndata:
    os.system('DatasetGen.py -d {} -s 200'.format(ns))
    t1.append(RunAverager(UseBlas=True))
    t2.append(RunAverager())


plt.figure()
plt.plot(ndata, t1, label='with blas')
plt.plot(ndata, t2, label='without blas')
plt.legend(numpoints=1, loc=0)
plt.xlabel('number of data points')
plt.ylabel('running time')
plt.savefig('TimevsData.pdf')