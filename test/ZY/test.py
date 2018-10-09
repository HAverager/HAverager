#!/usr/bin/env python

# Set path of the averager
import sys
sys.path.append('../bin')

# Import data reader and averager
import averager
import DataReader

# imports libraties for plotting
from numpy import *
import matplotlib.pyplot as plt
from matplotlib.cbook import get_sample_data
from operator import truediv

#read the data
data,stat,syst = DataReader.paverage('','','')

data = data*1000
stat = stat/100*data
syst = syst/100*data

#extract optional information
snames=DataReader.oerror
fnames=DataReader.fnames

#print input information
print 'Print Binning'
print DataReader.binnames

print 'Print Data'
print data

#print 'Print Syst'
#print syst

#print 'Print Stat'
#print stat

print snames

print 'Print Syst names'
trimsnames = map(str.strip, snames)
print trimsnames

#initialization (optional information)
averager.avin.initvariables()
averager.avin.setoutputfolder('./TOutP')
averager.avin.initeration = 0
averager.avin.setsnames(snames)
averager.avin.inwriteoriginal = True
averager.avin.infixstat = False
averager.avin.incorrectstatbias = False
averager.avin.inrescalestatsep = False

#perform averaging
dataAv,statAv,systAv = averager.average(data,stat,syst)

'''
#print output information
print "pull of data \n", averager.avout.pulldata

print "pull of systematics \n", averager.avout.pullsyst
print averager.avout.shiftsyst
print averager.avout.squeezesyst

print averager.avout.chi2
print averager.avout.ndof

#plotting

# Averaged data
# rotate initial arrays of data and stat uncertainty
swapdata=swapaxes(data,0,1)
swapstat=swapaxes(stat,0,1)

bins=arange(len(dataAv))
plt.figure()
fig, axs = plt.subplots(2, sharex=True)
axs[0].errorbar(bins, dataAv, yerr=statAv,ecolor='black',marker='o',ls='', label='Average')
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
	axs[0].errorbar(trimbins, trimdata, yerr=trimstat,ecolor='black',marker='o',ls='', label=fnames[ki-1])
axs[0].legend(numpoints=1)
axs[0].set_ylabel('data')
axs[0].set_ylim(-200, 4000)

# plot ratio plot
fig.subplots_adjust(hspace=0)
axs[1].errorbar(bins, map(truediv,dataAv,dataAv), 
	yerr=map(truediv,statAv,dataAv), ecolor='black',marker='o',ls='', label='Average')
for ki in range(len(swapdata)):
	axs[1].errorbar(bins, map(truediv,swapdata[ki-1],dataAv), yerr=map(truediv,statAv,dataAv), ecolor='black',marker='o',ls='', label=fnames[ki-1])

axs[1].set_ylim(0.6,1.4)
axs[1].set_xlim(-1,len(dataAv))
axs[1].set_xlabel('data points')
axs[1].set_ylabel('data / average')
#axs[1].legend(numpoints=1)

plt.plot()
plt.savefig('AvData.pdf')
plt.close()


# plot pulls of systematics
plt.figure()
plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.2)
plt.errorbar(arange(len(averager.avout.pullsyst)), 
	averager.avout.pullsyst, yerr=0,ecolor='black',marker='o',ls='')
plt.xlim(-1,len(averager.avout.pullsyst))
plt.ylabel('pull of systematics')
plt.xticks(arange(len(averager.avout.pullsyst)),trimsnames, rotation=45)
plt.setp(trimsnames)

plt.plot()
plt.savefig('AvPull.pdf')
plt.close()

# hist pulls of systematics
plt.figure()
plt.hist(averager.avout.pullsyst)
plt.ylabel("Systematic pulls")
plt.plot()
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


'''




