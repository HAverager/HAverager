#!/usr/bin/env python

import numpy as np
from optparse import OptionParser

# Dataset generator is a script to produce fake datasets in order to test averager

########################
# setup input parameters
########################
usage = ("DatasetGen.py parameters \nCreate dataset with gaussian uncertainties.  \n")
parser = OptionParser(usage)
parser.add_option("-d",dest="nData", default=50, type="int", help="Number of data points")
parser.add_option("-m",dest="nMes", default=2, type="int", help="Number of measurements")
parser.add_option("-s",dest="nSyst", default=10, type="int", help="Number of systematic sources")

parser.add_option("-M", dest="zMes", default=0.99, type="float", help="Fraction of non-emply measurements")
parser.add_option("-S", dest="zSyst", default=0.99, type="float", help="Fraction of non-emply systematics")

parser.add_option("--Poisson", action="store_true", dest="poisson", default=False, help="Use poissonian statistical uncertainties")
parser.add_option("--Smear", action="store_true", dest="smear", default=False, help="Systematic uncertainties are affected by statistics")

options, arguments = parser.parse_args()

# Number of data points, measurements, systematic
nData = options.nData
nMes = options.nMes
nSyst = options.nSyst

# part of non-emply measurement/systematics (between 0 and 1)
zMes = options.zMes
zSyst = options.zSyst

# Flags
doPoisson = options.poisson
doSmear = options.smear

# Parameters of uncertainties
# Minumal and maximal relative statistical/systematic uncertainty
vStat = [0.022,0.045]
vSyst = [0.03,0.08]



# array of truth data points
Tdata = nData*np.random.random_sample((nData))+1

# array of truth nuisanse parameters
Tshift = np.random.normal(0, 1, nSyst)

# array of statistical uncertainties. 
# 2D data points vs measurements
Mstat = (vStat[1]-vStat[0]) * np.random.random_sample((nMes, nData)) + vStat[0]

# array of systematic uncertainties. 
# 3D data points vs measurement vs systematic
# Systematics are linear across bins: Ax+B
base = np.ones((nSyst*nMes,nData))

Arnd = ((vSyst[1]-vSyst[0]) * np.random.random_sample(nSyst*nMes) + vSyst[0]).reshape(nSyst*nMes,1)
Brnd = (((vSyst[1]-vSyst[0]) * np.random.random_sample(nSyst*nMes) + vSyst[0])*np.sign(np.random.random_sample(nSyst*nMes)-0.5)).reshape(nSyst*nMes,1)

F = np.linspace(0.1,1.1,num=nData)

AXpB =  Arnd*base*F + Brnd
SystBase = AXpB.reshape(nMes,nSyst,nData)

if(doSmear):
	# Gaussian statistical fluctuation of the systematic uncertainties
	print 'Systematic uncertainties are smeared by statistics'
	Gaus3D  =  np.random.normal(0, 1, nMes*nData*nSyst).reshape((nMes, nSyst, nData))
	Msyst = (SystBase + ( Mstat.reshape(nMes,1,nData) * Gaus3D)).swapaxes(1,2)
else:
	# No statistical component in the systematics
	print 'Systematic uncertainties are not statistically fluctuated'
	Msyst = SystBase.swapaxes(1,2)

#Msyst = (vSyst[1]-vSyst[0]) * np.random.random_sample((nMes, nData, nSyst)) + vSyst[0]

# Add holes in the array of measurements and systematics
# - Some data points does not exists for a certain measurement
# - Some systematics does not exist for a certain measurement
# --------------------------------
Hdata = np.signbit(np.random.random_sample((nMes, nData))-zMes)*1
Hsyst = np.signbit(np.random.random_sample((nMes, 1, nSyst))-zSyst)*1

# array of measured data points
if (doPoisson):
	print 'Use Poissonian statistics \n'
	# Poisson statistics

	# size of statistics
	Vstat = 1/Mstat**2

	# Puisson fluctiation
	Pstat = np.random.poisson(Vstat)

	Mdata = Tdata*(1 + ((Pstat-Vstat)/Vstat) +np.sum(Msyst*Tshift*Hsyst,axis=2))*Hdata
else:
	print 'Use Gaussian statistics \n'
	# Gaussian statistics
	Gaus_ij = np.random.normal(0, 1, nMes*nData).reshape((nMes, nData))
	Mdata = Tdata*(1+(Mstat*Gaus_ij)+np.sum(Msyst*Tshift*Hsyst,axis=2))*Hdata

np.savetxt('Tshift.out', Tshift, fmt='%1.3f')
np.savetxt('Tdata.out', Tdata, fmt='%1.3f')

# Write output files. Python format
# Loop over measurements
for m in range(nMes):
	f = open('test'+str(m)+'.csv','w')
	f.write('bin1,data,stat')
	# Loop over systematics
	for s in range(nSyst):
		if(Hsyst[m][0][s]!=0):
			f.write(',error%05i' % s)
	f.write('\n')
	# Loop over data point
	for d in range(nData):
		if(np.abs(Mdata[m][d])>0.1):
			f.write('%4.0f,'% d)
			f.write('%8.3f,'% (Mdata[m][d]))
			f.write('%8.3f'% (Mstat[m][d]))#*Mdata[m][d]))	
			# Loop over systematics
			for s in range(nSyst):
				if(Hsyst[m][0][s]!=0):
					f.write(',%8.3f'% (Msyst[m][d][s]))#*Mdata[m][d]))
			f.write('\n')	
	f.close()

# Write output files. Fortran format
# Loop over measurements
for m in range(nMes):
	f = open('test'+str(m)+'.dat','w')

	f.write('&Data\n')
	f.write('   Name = \'Data%i\'\n'%nData)
	f.write('   NData = %i\n'% nData)
	f.write('   NColumn = %i\n'%(nSyst+3))
	f.write('   ColumnType = \'Bin\', \'Sigma\', %i*\'Error\'\n'%(nSyst+1))
	f.write('   ColumnName = \'Y\', \'x-section\', \'stat\'')
	for s in range(nSyst):
		if(Msyst[m][0][s]!=0):
			f.write(',\'error%05i\'' % s)
	f.write('\n')
	f.write('   Reaction = \'Bla\'\n')
	f.write('   Percent = false')
	for s in range(nSyst):
		if(Hsyst[m][0][s]!=0):
			f.write(',false')
	f.write('\n')
	f.write('&END\n')

	# Loop over data point
	for d in range(nData):
		if(Mdata[m][d]!=0):
			f.write('%4.0f '% d)
			f.write('%8.3f '% Mdata[m][d])		
			f.write('%8.3f '% (Mstat[m][d]*Mdata[m][d]))	
			# Loop over systematics
			for s in range(nSyst):
				if(Hsyst[m][0][s]!=0):
					f.write('%8.3f '% (Msyst[m][d][s]*Mdata[m][d]))
			f.write('\n')
	f.write('\n')
	f.close()

# Study correlation model
systBla = Msyst.reshape(nData*nMes,nSyst)

import matplotlib.pyplot as plt

Cov = (systBla).dot(systBla.transpose())
B  = np.diag(Cov)
nB = B.size
C = np.ones((nB,nB))
D = B.reshape(nB,1)*C*B
Corr = Cov/np.sqrt(D)


def PlotMatrix(matrix, fname, vmin=-1, vmax=1):
	plt.figure()
	im = plt.imshow(matrix, interpolation='none', alpha=None, vmin=vmin, vmax=vmax)
	plt.xlabel('bin number')
	plt.ylabel('bin number')
	clb = plt.colorbar()
	clb.set_label('Correlation', labelpad=-40, y=1.05, rotation=0)
	plt.savefig(fname)

PlotMatrix(Corr, fname='CorrGen.pdf')

