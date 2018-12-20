#!/usr/bin/env python

# Import data reader and averager
import PyAverager
import DataReader
import AveUtils
import AvePlot

import numpy as np

# read the data
data,stat,syst = DataReader.paverage('','','')

# Conver from percentrage to absolute
stat = stat/100*data
syst = syst/100*data

# extract optional information
snames=DataReader.oerror
fnames=DataReader.fnames

bins = DataReader.binnames
binnames = DataReader.bins

# Set input parameters
PyAverager.OutFolder = './TOutP'
PyAverager.nIterations = 0
PyAverager.writeoriginal = True
PyAverager.fixstat = False
PyAverager.correctstatbias = False
PyAverager.rescalestatsep = False
PyAverager.postrotatesyst = False
PyAverager.dosystimpact = False
PyAverager.ntoymc = 0
PyAverager.useblas = False

# Perform averaging
PyAverager.average(data,stat,syst,snames,fnames,bins,binnames)

# Compare results with fortran averager
Fout = 'TOutF'

# Differences in correlated matrises
CorrF = AveUtils.GetCorrMatrix(np.loadtxt(Fout+'/tab.dat')[:,1:])
CorrP = AveUtils.GetCorrMatrix(np.loadtxt(PyAverager.OutFolder+'/tab.dat')[:,1:])
AvePlot.PlotMatrix(CorrF-CorrP, 'TOutF.pdf', vmin=-0.01, vmax=0.01)

