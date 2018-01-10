#!/usr/bin/env python

# Set path of the averager
import sys
sys.path.append('../bin')
sys.path.append('../scripts')
# Import data reader and averager
import averager
import DataReader
import CorrModel

# imports libraties for plotting
from numpy import *
import matplotlib.pyplot as plt
from matplotlib.cbook import get_sample_data
from operator import truediv
from matplotlib.ticker import NullFormatter

def GetCorr(Cov):
	B  = diag(Cov)
	nB = B.size
	C = ones((nB,nB))
	D = B.reshape(nB,1)*C*B
	return Cov/sqrt(D)

def GetCov(A, axis=0):
	if(axis==0):
		Cov = A.dot(A.transpose())
	else:
		Cov = (A.transpose()).dot(A)
	return Cov


def PlotMatrix(matrix, fname, vmin=-1, vmax=1, title=''):
	plt.figure()
	im = plt.imshow(matrix, interpolation='none', alpha=None, vmin=vmin, vmax=vmax)
	plt.xlabel('bin number')
	plt.ylabel('bin number')
	plt.title(title)
	clb = plt.colorbar()
	clb.set_label('Correlation', labelpad=-40, y=1.05, rotation=0)
	plt.savefig(fname)

#perform averaging
def RunAverager(data,stat,syst, nToyMC=0, SysImp = False, itr=0, fixStat=False, corrStat=False):
	
	#initialization (optional information)
	averager.avin.cleaninvars()

	averager.avin.initvariables()
	averager.avin.setoutputfolder('./TOutP')
	averager.avin.initeration = itr
	#averager.avin.setsnames(snames)
	averager.avin.inwriteoriginal = True

	averager.avin.infixstat = False #fixStat
	averager.avin.incorrectstatbias = corrStat
	averager.avin.inrescalestatsep = fixStat#False

	averager.avin.inpostrotatesyst = False
	averager.avin.indosystimpact = SysImp
	averager.avin.inntoymc = nToyMC
	averager.avin.inuseblas = False

	dataAv,statAv,systAv = averager.average(data,stat,syst)
	return dataAv,statAv,systAv

#read the data
data,stat,syst = DataReader.paverage('','','')

syst=syst*data
stat=stat*data

# Printour ini Value
nBins = syst.shape[1]
nSyst = syst.shape[0]
nMes = syst.shape[2]


# Printour data for Blue
# allSyst = sqrt((syst**2).sum(axis=0) )
systBla = syst.transpose().reshape(nMes*nBins,nSyst)

PlotMatrix(GetCorr(GetCov(systBla)), fname='CorrIni.pdf', title='Input correlation model')

#extract optional information
snames=DataReader.oerror
fnames=DataReader.fnames

# Run default averager
dataAv,statAv,systAv = RunAverager(data,stat,syst, 0, False)
dataAv3,statAv3,systAv3 = RunAverager(data,stat,syst, 0, False, itr=5, fixStat=True)
dataAv4,statAv4,systAv4 = RunAverager(data,stat,syst, 0, False, itr=5, fixStat=True, corrStat=True)
dataAv5,statAv5,systAv5 = RunAverager(data,stat,syst, 0, False, itr=5, fixStat=False, corrStat=True)
dataAv2,statAv2,systAv2 = RunAverager(data,stat,syst, 0, False, itr=5)
totalAv = sqrt((systAv**2).sum(axis=0) + statAv**2 )

# Read Blue results
dataTrue = loadtxt('Tdata.out')
shiftTrue = loadtxt('Tshift.out')


# Plotting averaging result

# definitions for the axes
left, width = 0.1, 0.65
bottom, height = 0.1, 0.85
bottom_h = left_h = left + width + 0.01

rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, 0.2]
rect_histy = [left_h, bottom, 0.2, height]

# hist pulls of data
plt.figure()
axScatter = plt.axes(rect_scatter)
axScatter.errorbar(arange(len(averager.avout.pulldata[0])), 
	averager.avout.pulldata[0], yerr=0,ecolor='black',marker='o',ls='')
axScatter.set_xlim(-0.5,len(averager.avout.pulldata[0]))
plt.ylabel('pull of data')
plt.xlabel('index of measured point')
#plt.setp(trimsnames)
nullfmt = NullFormatter()
axHisty = plt.axes(rect_histy)
axHisty.yaxis.set_major_formatter(nullfmt)
axHisty.hist(averager.avout.pulldata[0], orientation='horizontal')
plt.xlabel('entry/bin')
plt.plot()
plt.savefig('AvPullData.pdf')
plt.close()


# plot pulls of systematics
plt.figure()
axScatter = plt.axes(rect_scatter)
axScatter.errorbar(arange(len(averager.avout.shiftsyst)), 
	nan_to_num(averager.avout.shiftsyst)*(-1), yerr=0,ecolor='black',marker='o',ls='')

axScatter.errorbar(arange(len(averager.avout.shiftsyst)), 
	shiftTrue, yerr=0,ecolor='black',marker='o',ls='')

plt.xlim(-0.5,len(averager.avout.shiftsyst))
#plt.ylim(-2,2)
plt.ylabel('pull of systematics')
plt.xlabel('index of systematic source')

axHisty = plt.axes(rect_histy)
axHisty.yaxis.set_major_formatter(nullfmt)
axHisty.hist(nan_to_num(averager.avout.pullsyst), orientation='horizontal')
plt.xlabel('entry/bin')
plt.plot()
plt.savefig('AvPullSyst.pdf')
plt.close()

# Averaged data
# rotate initial arrays of data and stat uncertainty
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


bins=arange(len(dataAv))
plt.figure()

plt.fill_between(bins, 1+(totalAv/dataTrue), 1-(totalAv/dataTrue),
    alpha=0.2, edgecolor='#1B2AFF', facecolor='#089FCC' )
plt.fill_between(bins, 1+(statAv/dataTrue), 1-(statAv/dataTrue),
    alpha=0.2, edgecolor='#1B2ACC', facecolor='#089FFF' )

plt.errorbar(bins, map(truediv,dataAv2,dataTrue), yerr=map(truediv,statAv,dataTrue), ecolor='black',marker='o',ls='', markersize='9', label='Mult corr')
plt.errorbar(bins, map(truediv,dataAv,dataTrue), yerr=map(truediv,statAv,dataTrue), ecolor='black',marker='o',ls='', label='Default')

plt.errorbar(bins, map(truediv,dataAv5,dataTrue), yerr=map(truediv,statAv,dataTrue), ecolor='black',marker='o',ls='', markersize='8', label='Stat corr')

plt.errorbar(bins, map(truediv,dataAv4,dataTrue), yerr=map(truediv,statAv,dataTrue), ecolor='black',marker='o',ls='', label='Stat corr + FixStat')

plt.errorbar(bins, map(truediv,swapdata[0],dataTrue), yerr=map(truediv,statAv,dataTrue), ecolor='black',marker='o',ls='', label='Dataset 1')

plt.errorbar(bins, map(truediv,swapdata[1],dataTrue), yerr=map(truediv,statAv,dataTrue), ecolor='black',marker='o',ls='', label='Dataset 2')


plt.legend(numpoints=1, loc=0)
plt.ylim(0.1,1.9)
plt.xlabel('data points')
plt.ylabel('All / Data True')
plt.plot()
plt.savefig('AvData.pdf')
plt.close()


'''
Closure - Smear Stat, Gauss
Closure 2 - Gaus, No Smear
3 - Poisson, Smear
4 Poisson, No Smear

'''

