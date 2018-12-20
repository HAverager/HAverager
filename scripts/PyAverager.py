#!/usr/bin/env python

import averager

import numpy as np
import pandas as pd
import os
import glob

# Set default input parameters
OutFolder = './TOutP'
nIterations = 0
writeoriginal = True
fixstat = False
correctstatbias = False
rescalestatsep = False
postrotatesyst = False
dosystimpact = False
ntoymc = 0
useblas = False
debug = 0
nData = 0


def getDummyNames(nnames, templ):
    names = np.chararray(nnames, itemsize=15)
    names[:] = templ
    names = names + np.arange(nnames).astype('str')
    names = [x.ljust(32) for x in names]
    #names = list(np.core.defchararray.ljust(names, 32))
    return names


def average(data, stat, syst, snames=None, fnames=None, bins=None, binnames=None):
    # initialization
    averager.avin.initvariables()
    averager.avin.setoutputfolder(OutFolder)
    averager.avin.initeration = nIterations
    averager.avin.inwriteoriginal = writeoriginal
    averager.avin.indebug = debug

    averager.avin.infixstat = fixstat
    averager.avin.incorrectstatbias = correctstatbias
    averager.avin.inrescalestatsep = rescalestatsep

    averager.avin.inpostrotatesyst = postrotatesyst
    averager.avin.indosystimpact = dosystimpact
    averager.avin.inntoymc = ntoymc
    averager.avin.inuseblas = useblas

    # set internal averager parameters
    if bins is None:
        bins = np.arange(data.shape[0]+1)

    if binnames is None:
        ndim = 1 if len(bins.shape)<2 else bins.shape[1]
        binnames = getDummyNames(ndim, 'Bin')
    else:
        binnames = [x.ljust(32) for x in binnames]

    averager.avin.setbinning(bins, binnames)

    if snames is None:
        snames = getDummyNames(syst.shape[0], 'syst')
    else:
        snames = [x.ljust(32) for x in snames]

    if fnames is None:
        fnames = getDummyNames(data.shape[1], 'file')
    else:
        fnames = [x.ljust(32) for x in fnames]

    # Perform averaging
    dataAv, statAv, systAv = averager.average(data, stat, syst, snames, fnames)

    # Define global output variables
    global nData
    global nSyst
    global nFiles
    global chi2
    global ndof

    nData = data.shape[0]
    nFiles = data.shape[1]
    nSyst = syst.shape[0]

    chi2, ndof = averager.getchi2()

    return dataAv, statAv, systAv

def getSystPulls():
    return averager.getpulls(nSyst)

def getDataPulls():
    return averager.getdatapulls(nData, nFiles)

def getSystImpact():
    return averager.getsysimpact(nData, nSyst)

def getSystShifts():
    return averager.getshiftsyst(nSyst)

def getToyStat():
    return averager.gettoystat(nData)
