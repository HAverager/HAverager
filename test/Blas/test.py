#!/usr/bin/env python

import os
import time
import matplotlib.pyplot as plt

# Import data reader and averager
import PyAverager
import DataReader

# Script to test time performance of HAvergaer with different options

# perform averaging
def RunAverager(UseBlas=False, nToyMC=0, SysImp=False, itr=0):

    # read the data
    print 'read data'
    data, stat, syst = DataReader.paveragefast()
    print 'do Averaging'

    # Set input parameters
    PyAverager.OutFolder = './TOutP'
    PyAverager.nIterations = itr
    PyAverager.writeoriginal = False

    PyAverager.fixstat = False
    PyAverager.correctstatbias = False
    PyAverager.rescalestatsep = False

    PyAverager.postrotatesyst = False
    PyAverager.dosystimpact = SysImp
    PyAverager.ntoymc = nToyMC
    PyAverager.useblas = UseBlas

    # Perform averaging
    start = time.time()
    PyAverager.average(data, stat, syst)
    end = time.time()

    return (end-start)

# Test 1.
nsyst = [20, 50, 100, 150, 200, 250, 300, 450, 1000, 1500, 2000]
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
ndata = [100, 200, 500, 1000, 1400]
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
