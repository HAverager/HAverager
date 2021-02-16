#!/usr/bin/env python

# Import data reader and averager
import PyAverager
import DataReader
import AveUtils
import AvePlot

# imports libraties for plotting
from numpy import *

#perform averaging
def RunAverager(data,stat,syst, nToyMC=0, SysImp = False, itr=0, fixStat=False, corrStat=False):
	
	# Set input parameters
	PyAverager.OutFolder = './TOutP'
	PyAverager.nIterations = 0
	PyAverager.writeoriginal = True

	PyAverager.fixstat = False
	PyAverager.correctstatbias = corrStat
	PyAverager.rescalestatsep = fixStat

	PyAverager.postrotatesyst = False
	PyAverager.dosystimpact = SysImp
	PyAverager.ntoymc = nToyMC
	PyAverager.useblas = False

	# Perform averaging
	dataAv,statAv,systAv = PyAverager.average(data, stat, syst)
	return dataAv,statAv,systAv

# read the data
data,stat,syst = DataReader.paverage('','','')

nBins = syst.shape[1]
nSyst = syst.shape[0]
nMes = syst.shape[2]

# Plot correlation matrix
systBla = syst.transpose().reshape(nMes*nBins,nSyst)
AvePlot.PlotMatrix(AveUtils.GetCorrMatrix(systBla), fname='CorrIni.pdf', title='Input correlation model')

# Run averager with different options
dataAv,statAv,systAv = RunAverager(data,stat,syst, 0, False)
dataAv3,statAv3,systAv3 = RunAverager(data,stat,syst, 0, False, itr=5, fixStat=True)
dataAv4,statAv4,systAv4 = RunAverager(data,stat,syst, 0, False, itr=5, fixStat=True, corrStat=True)
dataAv5,statAv5,systAv5 = RunAverager(data,stat,syst, 0, False, itr=5, fixStat=False, corrStat=True)
dataAv2,statAv2,systAv2 = RunAverager(data,stat,syst, 0, False, itr=5)
totalAv = sqrt((systAv**2).sum(axis=0) + statAv**2 )

# Read truth information
dataTrue = loadtxt('Tdata.out')
shiftTrue = loadtxt('Tshift.out')


# Plotting averaging result
AvePlot.PlotHistDecompose(PyAverager.getDataPulls()[0], 'AvPullData.pdf',
						  ytitle='pull of data', xtitle='index of measured point')

AvePlot.PlotHistDecompose(PyAverager.getSystShifts(), 'AvPullSyst.pdf',
				  ytitle='pull of systematics', xtitle='index of systematic source')



# Plot averaged data
swapdata=swapaxes(data,0,1)
swapstat=swapaxes(stat,0,1)

# Use last bin as averge over bins
dataTrue[len(dataTrue)-1] = dataTrue.mean()
dataAv[len(dataAv)-1] = dataAv.mean()
dataAv2[len(dataAv2)-1] = dataAv2.mean()
dataAv3[len(dataAv3)-1] = dataAv3.mean()
dataAv4[len(dataAv4)-1] = dataAv4.mean()
dataAv5[len(dataAv5)-1] = dataAv5.mean()
swapdata[0][len(swapdata[0])-1] = swapdata[0].mean()
swapdata[1][len(swapdata[1])-1] = swapdata[1].mean()
statAv[len(dataTrue)-1] = statAv[len(dataTrue)-1]/sqrt(len(dataTrue)-1)

dataAll = stack((dataTrue, dataAv2, dataAv, dataAv5, dataAv4, swapdata[0], swapdata[1]))
labels = ['Data True', 'Mult corr', 'Default', 'Stat corr', 'Stat corr + FixStat', 'Dataset 1', 'Dataset 2']

AvePlot.PlotProfile(dataAll, statAv, totalAv, labels, 'AvData.pdf',
			ratio=2, xtitle='data points', ytitle='All / Data True', rlims=[0.1, 1.9])
