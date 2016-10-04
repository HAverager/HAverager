#!/usr/bin/env python

import numpy as np

# Dataset generator is a script to produce fake datasets in order to test averager

########################
# setup input parameters
########################

# Number of data points, measurements, systematic
nData = 10
nMes = 2
nSyst = 5

# part of non-emply measurement/systematics (between 0 and 1)
zMes = 0.99 
zSyst = 0.99

# Min and mix number of measurements for data point (have to be smaller as total number of measurements) 
#rMes = [1,3]

# Min and mix number of systematic uncertainties for measurement (have to be smaller as total number of systematic uncertainties) 
#rSyst = [3,4]


# Parameters of uncertainties
# Minumal and maximal relative statistical/systematic uncertainty
vStat = [0.01,0.03]
vSyst = [0.08,0.11]



# array of truth data points
Tdata = nData*np.random.random_sample((nData))

# array of truth nuisanse parameters
Tshift = np.random.normal(0, 1, nSyst)

# array of statistical uncertainties. 
# 2D data points vs measurements
Mstat = (vStat[1]-vStat[0]) * np.random.random_sample((nMes, nData)) + vStat[0]

# array of systematic uncertainties. 
# 3D data points vs measurement vs systematic
Msyst = (vSyst[1]-vSyst[0]) * np.random.random_sample((nMes, nData, nSyst)) + vSyst[0]

# Add holes in the array of measurements and systematics
# - Some data points does not exists for a certain measurement
# - Some systematics does not exist for a certain measurement
# --------------------------------
Hdata = np.signbit(np.random.random_sample((nMes, nData))-zMes)*1
Hsyst = np.signbit(np.random.random_sample((nMes, 1, nSyst))-zSyst)*1

# array of measured data points
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
		if(Mdata[m][d]!=0):
			f.write('%4.0f,'% d)
			f.write('%8.3f,'% (Mdata[m][d]))		
			f.write('%8.3f'% (Mstat[m][d]*Mdata[m][d]))	
			# Loop over systematics
			for s in range(nSyst):
				if(Hsyst[m][0][s]!=0):
					f.write(',%8.3f'% (Msyst[m][d][s]*Mdata[m][d]))
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
