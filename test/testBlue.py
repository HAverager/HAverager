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


def PlotMatrix(matrix, fname, vmin=-1, vmax=1):
	plt.figure()
	im = plt.imshow(matrix, interpolation='none', alpha=None, vmin=vmin, vmax=vmax)
	plt.xlabel('bin number')
	plt.ylabel('bin number')
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

# Printour data for Blue
allSyst = sqrt((syst**2).sum(axis=0) )
systBla = syst.transpose().reshape(43*2,54)

savetxt('ZeePtBlue2.txt', column_stack((data[:,0],stat[:,0],allSyst[:,0])))
savetxt('ZmmPtBlue2.txt', column_stack((data[:,1],stat[:,1],allSyst[:,1])))

PlotMatrix(corrcoef( systBla.dot(systBla.transpose()) ), fname='CorrIni.pdf')
savetxt('ZPtBlueCor.txt', corrcoef( systBla.dot(systBla.transpose()) ), fmt='%.7e',)

#extract optional information
snames=DataReader.oerror
fnames=DataReader.fnames

# Run default averager
dataAv,statAv,systAv = RunAverager(data,stat,syst, 200, True)

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

BlueCov = BlueCov + diag(statBlue**2)

print 'Data Diff \n', (dataAv-dataBlue)/dataAv*100
print 'Total Diff \n',(totalAv-totalBlue)/totalAv*100
print 'TotalToy Diff \n',(totaltoy-totalBlue)/totaltoy*100
print 'Stat Diff \n',(Toystat-statBlue)/Toystat*100
print 'Syst Diff \n',(totalSystImpact-systBlue)/totalSystImpact*100

PlotMatrix(corrcoef(CovAv), fname='CorrHave.pdf')
PlotMatrix(corrcoef(BlueCov), fname='CorrBlue.pdf')
PlotMatrix(corrcoef(BlueCov)-corrcoef(CovAv), fname='CorrDiffBlue.pdf')
PlotMatrix(corrcoef(CovAvImpact)-corrcoef(CovAv), fname='CorrDiffImpact.pdf')

# Reduce correlation model of input information
syst2, stat2 = CorrModel.ReduceCorrelationModel(syst, coef=0.01)
stat2 = sqrt(stat**2 + stat2**2)


# perform averaging
dataAv2,statAv2,systAv2 = RunAverager(data,stat2,syst2, 0, False)
totalAv2 = sqrt((systAv2**2).sum(axis=0) + statAv2**2)
CovAv2 = (systAv2.transpose()).dot(systAv2) + diag(statAv2**2)

print (dataAv2-dataAv)/statAv
print (statAv2-statAv)/statAv2
print (totalAv2-totalAv)/totalAv2


# Further reduction of correlation model of input information
syst3, stat3 = CorrModel.ReduceCorrelationModel(syst, stat, coef=0.01)

# Perform averaging
dataAv3,statAv3,systAv3 = RunAverager(data,stat3,syst3, 0, False)
totalAv3 = sqrt((systAv3**2).sum(axis=0) + statAv3**2)
CovAv3 = (systAv3.transpose()).dot(systAv3) + diag(statAv3**2)

print (totalAv3-totalAv)/totalAv3

# Reduce correlation for output
systAv4, statAv4 = CorrModel.ReduceCorrelationModel(systAv3, coef=0.4)
CovAv4 = (systAv4.transpose()).dot(systAv4) + diag(statAv4**2) + diag(statAv3**2)

PlotMatrix(corrcoef(CovAv)-corrcoef(CovAv2), fname='CorrInput.pdf')
PlotMatrix(corrcoef(CovAv)-corrcoef(CovAv3), fname='CorrInputStat.pdf')
PlotMatrix(corrcoef(CovAv3)-corrcoef(CovAv4), fname='CorrOutput.pdf')


# Plot results
# Total uncertainty
plt.figure()
plt.plot((totalAv3-totalAv)/totalAv, label='Input stat')
plt.plot((totalAv2-totalAv)/totalAv, label='Input')
plt.plot((totalBlue-totalAv)/totalAv, label='Blue')
plt.plot((totaltoy-totalAv)/totalAv, label='Toy')
plt.legend(loc=2)
plt.xlabel('bin number')
plt.ylabel('(X-Total) / Total')
#plt.yscale('log')
plt.savefig('TotalUncertainty.pdf')


plt.savefig('AvPullHist.pdf')
plt.close()

# hist pulls of data
for i in range(len(averager.avout.pulldata)):
	plt.figure()
	plt.hist(averager.avout.pulldata[i])
	plt.ylabel("Data pulls")
	plt.plot()

	plt.savefig('AvDataPullHist'+str(i+1)+'.pdf')
	plt.close()

# hist pulls of data
plt.figure()
plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.2)
plt.errorbar(arange(len(averager.avout.pulldata[0])), 
	averager.avout.pulldata[0], yerr=0,ecolor='black',marker='o',ls='')
plt.xlim(-1,len(averager.avout.pulldata[0]))
plt.ylabel('pull of data')
plt.setp(trimsnames)

plt.plot()
plt.savefig('AvPullData.pdf')
plt.close()


