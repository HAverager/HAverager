#!/usr/bin/env python

# Import libraries
import numpy as np
import matplotlib.pyplot as plt
import sys

# read data from txt files
folder = sys.argv[1]+'/'

systpull = np.loadtxt(folder+'sys.txt', usecols=[4])
systnames = np.loadtxt(folder+'sys.txt', dtype='str', usecols=[1])
data = np.loadtxt(folder+'tab.dat')
ndata = len(data)
pulldata = np.loadtxt(folder+'chi2map.dat', usecols=[3])
#pulldata = pulldata[0:ndata]

print 'Number of data points: ', ndata
print 'Number of measurements: ', len(pulldata)/ndata 
print 'Number of systematics: ', len(systnames)
print 'List of systematics: ', systnames


# plotting

# hist pulls of systematics
plt.figure()
plt.hist(systpull)
plt.xlabel("Systematic pulls")
plt.plot()
plt.savefig('AvPullHist.pdf')
plt.close()


# hist pulls of data
for i in range(len(pulldata)/ndata):
	cutpulldata = pulldata[ndata*i:ndata*(i+1)]
	plt.figure()
	plt.hist(cutpulldata)
	plt.xlabel("Data pulls")
	plt.plot()
	plt.savefig('AvDataPullHist'+str(i+1)+'.pdf')
	plt.close()


