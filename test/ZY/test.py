#!/usr/bin/env python

# Import data reader and averager
import averager
import DataReader

#read the data
data,stat,syst = DataReader.paverage('','','')

print stat
stat = stat/100*data
syst = syst/100*data

# extract optional information
snames=[x.ljust(32) for x in DataReader.oerror]
fnames=[x.ljust(32) for x in DataReader.fnames]
bins = [x.ljust(32) for x in DataReader.bins]
binnames = DataReader.binnames

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

#initialization (optional information)
averager.avin.initvariables()
averager.avin.setoutputfolder('./TOutP')
averager.avin.initeration = 0
averager.avin.indebug = 2
averager.avin.setsnames(snames)
averager.avin.setfnames(fnames)
averager.avin.setbinning(binnames, bins)
averager.avin.inwriteoriginal = True
averager.avin.infixstat = False
averager.avin.incorrectstatbias = False
averager.avin.inrescalestatsep = True

#perform averaging
dataAv,statAv,systAv = averager.average(data,stat,syst)