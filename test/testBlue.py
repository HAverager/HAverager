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
def RunAverager(data,stat,syst, nToyMC=0, SysImp = False):
	
	#initialization (optional information)
	averager.avin.cleaninvars()

	averager.avin.initvariables()
	averager.avin.setoutputfolder('./TOutP')
	averager.avin.initeration = 0
	#averager.avin.setsnames(snames)
	averager.avin.inwriteoriginal = True
	averager.avin.infixstat = False
	averager.avin.incorrectstatbias = False
	averager.avin.inrescalestatsep = False
	averager.avin.inpostrotatesyst = False

	averager.avin.indosystimpact = SysImp
	averager.avin.inntoymc = nToyMC
	averager.avin.inuseblas = False

	dataAv,statAv,systAv = averager.average(data,stat,syst)
	return dataAv,statAv,systAv

#read the data
data,stat,syst = DataReader.paverage('','','')

# Scale data
data = data*1000
stat = stat/100*data
syst = syst/100*data

print syst.shape
nBins = syst.shape[1]
nSyst = syst.shape[0]

# Printour data for Blue
allSyst = sqrt((syst**2).sum(axis=0) )
systBla = syst.transpose().reshape(43*2,nSyst)

savetxt('ZeePtBlue2.txt', column_stack((data[:,0],stat[:,0],allSyst[:,0])))
savetxt('ZmmPtBlue2.txt', column_stack((data[:,1],stat[:,1],allSyst[:,1])))

PlotMatrix(GetCorr(GetCov(systBla)), fname='CorrIni.pdf', title='Input correlation model')
savetxt('ZPtBlueCor.txt', GetCorr(GetCov(systBla)), fmt='%.7e',)


#extract optional information
snames=DataReader.oerror
fnames=DataReader.fnames

# Run default averager
dataAv,statAv,systAv = RunAverager(data,stat,syst, 200, True)

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
	nan_to_num(averager.avout.pullsyst), yerr=0,ecolor='black',marker='o',ls='')
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

'''
# Averaged data
# rotate initial arrays of data and stat uncertainty
swapdata=swapaxes(data,0,1)
swapstat=swapaxes(stat,0,1)
hlabel = ['Dataset 1', 'Dataset 2']
bins=arange(len(dataAv))
plt.figure()
fig, axs = plt.subplots(2, sharex=True)

# Loop over all data samples
for ki in range(len(swapdata)):
#   Remove 0s from arrays
	trimbins = []
	trimdata = []
	trimstat = []
	for i in range(len(bins)):
		if(swapdata[ki-1][i]!=0.0):
			trimbins.append(bins[i])
			trimdata.append(swapdata[ki-1][i])
			trimstat.append(swapstat[ki-1][i])
	axs[0].errorbar(trimbins, trimdata, yerr=trimstat,ecolor='black',marker='o',ls='', label=hlabel[ki])
axs[0].errorbar(bins, dataAv, yerr=statAv,ecolor='black',marker='o',ls='', label='Average')
axs[0].legend(numpoints=1, loc=0)
axs[0].set_ylabel('data')
axs[0].set_ylim(1e-5, 1e+2)
axs[0].set_yscale('log')

# plot ratio plot
fig.subplots_adjust(hspace=0)
#axs[1].errorbar(bins, map(truediv,dataAv,dataAv), 
#	yerr=map(truediv,statAv,dataAv), ecolor='black',marker='o',ls='', label='Average')
axs[1].fill_between(bins, 1+(statAv/dataAv), 1-(statAv/dataAv),
    alpha=0.2, edgecolor='#1B2ACC', facecolor='#089FFF' )

for ki in range(len(swapdata)):
	axs[1].errorbar(bins, map(truediv,swapdata[ki-1],dataAv), yerr=map(truediv,statAv,dataAv), ecolor='black',marker='o',ls='', label=fnames[ki-1])

axs[1].set_ylim(0.94,1.06)
axs[1].set_xlim(-1,len(dataAv))
axs[1].set_xlabel('data points')
axs[1].set_ylabel('data / average')
#axs[1].legend(numpoints=1)

plt.plot()
plt.savefig('AvData.pdf')
plt.close()
'''

Toystat = averager.avout.toystat
SystImpact = averager.avout.sysimpact
totalAv = sqrt((systAv**2).sum(axis=0) + statAv**2 )
totalSystImpact = sqrt((SystImpact**2).sum(axis=0) )
totaltoy = sqrt(totalSystImpact**2 + Toystat**2)
CovAv = (systAv.transpose()).dot(systAv) + diag(statAv**2)
CovAvImpact = (SystImpact.transpose()).dot(SystImpact) + diag(Toystat**2)

# Comparison with Blue
# Read Blue results
BlueVals = loadtxt('ZptBlueOut.txt')
BlueCov = loadtxt('ZptBlueCov.txt')

dataBlue = BlueVals[:,0]
totalBlue = sqrt((BlueVals[:,1:]**2).sum(axis=1) )
statBlue = BlueVals[:,1]
systBlue = BlueVals[:,2]

BlueCov = BlueCov

print 'Data Diff \n', (dataAv-dataBlue)/dataAv*100
print 'Total Diff \n',(totalAv-totalBlue)/totalAv*100
print 'TotalToy Diff \n',(totaltoy-totalBlue)/totaltoy*100
print 'Stat Diff \n',(Toystat-statBlue)/Toystat*100
print 'Syst Diff \n',(totalSystImpact-systBlue)/totalSystImpact*100



# Reduce correlation model of input information
syst2, stat2 = CorrModel.ReduceCorrelationModel(syst, coef=0.01)
stat2 = sqrt(stat**2 + stat2**2)


# perform averaging
dataAv2,statAv2,systAv2 = RunAverager(data,stat2,syst2, 0, False)
totalAv2 = sqrt((systAv2**2).sum(axis=0) + statAv2**2)
CovAv2 = (systAv2.transpose()).dot(systAv2) + diag(statAv2**2)


# Further reduction of correlation model of input information
syst3, stat3 = CorrModel.ReduceCorrelationModel(syst, stat, coef=0.01)

# Perform averaging
dataAv3,statAv3,systAv3 = RunAverager(data,stat3,syst3, 0, False)
totalAv3 = sqrt((systAv3**2).sum(axis=0) + statAv3**2)
CovAv3 = (systAv3.transpose()).dot(systAv3) + diag(statAv3**2)


# Reduce correlation for output
systAv4, statAv4 = CorrModel.ReduceCorrelationModel(systAv3, coef=0.2)
CovAv4 = (systAv4.transpose()).dot(systAv4) + diag(statAv4**2) + diag(statAv3**2)


# Plot results
PlotMatrix(GetCorr(CovAv)-GetCorr(CovAv2), fname='CorrDiffInput.pdf', vmin=-0.2, vmax=0.2, title='Input syst. redused by 1%')
PlotMatrix(GetCorr(CovAv)-GetCorr(CovAv3), fname='CorrDiffInputStat.pdf', title='Input syst.+stat. reduset by 1%')
PlotMatrix(GetCorr(CovAv3)-GetCorr(CovAv4), fname='CorrDiffOutput.pdf', title='Output syst.+stat. reduced by 20%')
PlotMatrix(GetCorr(CovAv), fname='CorrHave.pdf')
PlotMatrix(GetCorr(BlueCov), fname='CorrBlue.pdf')
PlotMatrix(GetCorr(BlueCov)-GetCorr(CovAv), fname='CorrDiffBlue.pdf', vmin=-0.2, vmax=0.2, title='Blue-Averager(Default)')
PlotMatrix(GetCorr(CovAvImpact)-GetCorr(CovAv), fname='CorrDiffImpact.pdf', vmin=-0.2, vmax=0.2, title='Averager(Impact)-Averager(Default)')
PlotMatrix(GetCorr(BlueCov)-GetCorr(CovAvImpact), fname='CorrDiffImpactBlue.pdf', vmin=-0.2, vmax=0.2, title='Blue-Averager(Impact)')


# Total uncertainty
plt.figure()
plt.plot((totalAv3-totalAv)/totalAv, label='Input stat')
plt.plot((totalAv2-totalAv)/totalAv, label='Input')
plt.plot((totalBlue-totalAv)/totalAv, label='Blue')
plt.plot((totaltoy-totalAv)/totalAv, label='Impact')
plt.legend(loc=0)
plt.xlabel('bin number')
plt.ylabel('(X-Total) / Total')
#plt.yscale('log')
plt.savefig('TotalUncertainty.pdf')


# Central values
plt.figure()
plt.plot((dataAv3-dataAv)/dataAv, label='Input stat')
plt.plot((dataAv2-dataAv)/dataAv, label='Input')
plt.plot((dataBlue-dataAv)/dataAv, label='Blue')
#plt.plot((totaltoy-dataAv)/dataAv, label='Toy')
plt.legend(loc=0)
plt.xlabel('bin number')
plt.ylabel('(X-default) / default')
#plt.yscale('log')
plt.savefig('CentralValues.pdf')



