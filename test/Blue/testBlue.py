#!/usr/bin/env python

# Import data reader and averager
import PyAverager
import DataReader
import AveUtils
import AvePlot

# imports libraties for plotting
from numpy import *
import matplotlib.pyplot as plt


#perform averaging
def RunAverager(data,stat,syst, nToyMC=0, SysImp = False):

	# Set input parameters
	PyAverager.OutFolder = './TOutP'
	PyAverager.nIterations = 0
	PyAverager.writeoriginal = True
	PyAverager.fixstat = False
	PyAverager.correctstatbias = False
	PyAverager.rescalestatsep = False
	PyAverager.postrotatesyst = False
	PyAverager.dosystimpact = SysImp
	PyAverager.ntoymc = nToyMC
	PyAverager.useblas = False

	# Perform averaging
	dataAv,statAv,systAv = PyAverager.average(data, stat, syst)
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

AvePlot.PlotMatrix(AveUtils.GetCorrMatrix(systBla), fname='CorrIni.pdf', title='Input correlation model')
savetxt('ZPtBlueCor.txt', AveUtils.GetCorrMatrix(systBla), fmt='%.7e',)


#extract optional information
snames=DataReader.oerror
fnames=DataReader.fnames

# Run default averager
dataAv,statAv,systAv = RunAverager(data,stat,syst, 200, True)

Toystat = PyAverager.getToyStat()
SystImpact = PyAverager.getSystImpact()
totalAv = sqrt((systAv**2).sum(axis=0) + statAv**2 )
totalSystImpact = sqrt((SystImpact**2).sum(axis=0) )
totaltoy = sqrt(totalSystImpact**2 + Toystat**2)
CovAv = (systAv.transpose()).dot(systAv) + diag(statAv**2)
CovAvImpact = (SystImpact.transpose()).dot(SystImpact) + diag(Toystat**2)

# Plotting averaging result
AvePlot.PlotHistDecompose(PyAverager.getDataPulls()[0], 'AvDataPullHist2.pdf',
						  ytitle='pull of data', xtitle='index of measured point')

AvePlot.PlotHistDecompose(PyAverager.getSystShifts(), 'AvPullHist2.pdf',
				  ytitle='pull of systematics', xtitle='index of systematic source')


swapdata=swapaxes(data,0,1)
swapstat=swapaxes(stat,0,1)
hlabel = ['Averaged', 'Dataset 1', 'Dataset 2']
dataAll = stack((dataAv, swapdata[0], swapdata[1]))

AvePlot.PlotProfile(dataAll, statAv, totalAv, hlabel, 'AvData.pdf',
			ratio=1, xtitle='data points', ytitle='All / Averaged data', rlims=[0.94, 1.06])



# Comparison with Blue
# Read Blue results
BlueVals = loadtxt('ZptBlueOut.txt')
BlueCov = loadtxt('ZptBlueCov.txt')

dataBlue = BlueVals[:,0]
totalBlue = sqrt((BlueVals[:,1:]**2).sum(axis=1) )
statBlue = BlueVals[:,1]
systBlue = BlueVals[:,2]

# Reduce correlation model of input information
syst2, stat2 = AveUtils.ReduceCorrelationModel(syst, coef=0.01)
stat2 = sqrt(stat**2 + stat2**2)

# perform averaging
dataAv2,statAv2,systAv2 = RunAverager(data,stat2,syst2, 0, False)
totalAv2 = sqrt((systAv2**2).sum(axis=0) + statAv2**2)
CovAv2 = (systAv2.transpose()).dot(systAv2) + diag(statAv2**2)

# Further reduction of correlation model of input information
syst3, stat3 = AveUtils.ReduceCorrelationModel(syst, stat, coef=0.01)

# Perform averaging
dataAv3,statAv3,systAv3 = RunAverager(data,stat3,syst3, 0, False)
totalAv3 = sqrt((systAv3**2).sum(axis=0) + statAv3**2)
CovAv3 = (systAv3.transpose()).dot(systAv3) + diag(statAv3**2)

# Reduce correlation for output
systAv4, statAv4 = AveUtils.ReduceCorrelationModel(systAv3, coef=0.2)
CovAv4 = (systAv4.transpose()).dot(systAv4) + diag(statAv4**2) + diag(statAv3**2)

# Plot results
AvePlot.PlotMatrix(AveUtils.GetCorr(CovAv)-AveUtils.GetCorr(CovAv2), fname='CorrDiffInput.pdf', vmin=-0.2, vmax=0.2, title='Input syst. redused by 1%')
AvePlot.PlotMatrix(AveUtils.GetCorr(CovAv)-AveUtils.GetCorr(CovAv3), fname='CorrDiffInputStat.pdf', title='Input syst.+stat. reduset by 1%')
AvePlot.PlotMatrix(AveUtils.GetCorr(CovAv3)-AveUtils.GetCorr(CovAv4), fname='CorrDiffOutput.pdf', title='Output syst.+stat. reduced by 20%')
AvePlot.PlotMatrix(AveUtils.GetCorr(CovAv), fname='CorrHave.pdf')
AvePlot.PlotMatrix(AveUtils.GetCorr(BlueCov), fname='CorrBlue.pdf')
AvePlot.PlotMatrix(AveUtils.GetCorr(BlueCov)-AveUtils.GetCorr(CovAv), fname='CorrDiffBlue.pdf', vmin=-0.2, vmax=0.2, title='Blue-Averager(Default)')
AvePlot.PlotMatrix(AveUtils.GetCorr(CovAvImpact)-AveUtils.GetCorr(CovAv), fname='CorrDiffImpact.pdf', vmin=-0.2, vmax=0.2, title='Averager(Impact)-Averager(Default)')
AvePlot.PlotMatrix(AveUtils.GetCorr(BlueCov)-AveUtils.GetCorr(CovAvImpact), fname='CorrDiffImpactBlue.pdf', vmin=-0.2, vmax=0.2, title='Blue-Averager(Impact)')


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

