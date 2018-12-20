#!/usr/bin/env python

import numpy as np
from optparse import OptionParser
import AvePlot
import AveUtils

# Write output files. CSV format
def WriteCSV(data, syst, stat, doPercent = False):

	# Scale to correct representation
	if doPercent:
		stat = stat*100
		syst = syst*100
	else:
		stat = stat*data
		syst = syst*data[:,:,None]

	# Create array of all data
	dataAll = np.concatenate((data[:,:,None], stat[:,:,None], syst ), axis=2)

	# Create binning
	bins = np.arange(data.shape[1])[None,:].repeat(data.shape[0], axis=0)

	# Create array of all data
	dataAll = np.concatenate((bins[:,:,None],
							  data[:,:,None], stat[:,:,None], syst ), axis=2)


	# Create header
	heds = np.chararray(syst.shape[2]+3, itemsize=15)
	heds[0] = 'bin1'
	heds[1] = 'data'
	heds[2] = 'stat'
	heds[3:] = 'error'
	heds[3:] = heds[3:] + np.arange(syst.shape[2]).astype('str')

	# Loop over files
	for i, fdata in enumerate(dataAll):

		# remove columns and rows with 0
		fdatar  = np.delete(fdata, np.where(~fdata[:,1:].any(axis=1)), axis=0)
		fdatacr = np.delete(fdatar, np.where(~fdatar.any(axis=0)), axis=1)

		# Create header
		header = (np.delete(heds, np.where(~fdatar.any(axis=0)), axis=0)).tolist()

		# Save file
		np.savetxt('test'+str(i)+'.csv', fdatacr, delimiter=',', fmt='%.5f',
			   header=','.join(header), comments='')

# Write output files. Fortran format
def writeDAT(Mdata, Msyst, Mstat, Hsyst, doPercent):

	nData = Mdata.shape[1]
	nFiles = Mdata.shape[0]
	nSyst = Msyst.shape[2]

	# Loop over measurements
	for m in range(nFiles):
		f = open('test' + str(m) + '.dat', 'w')

		f.write('&Data\n')
		f.write('   Name = \'Data%i\'\n' % nData)
		f.write('   NData = %i\n' % (Mdata[m, :] != 0).sum())
		f.write('   NColumn = %i\n' % (Hsyst[m, 0, :].sum() + 3))
		f.write('   ColumnType = \'Bin\', \'Sigma\', %i*\'Error\'\n' % (Hsyst[m, 0, :].sum() + 1))
		f.write('   ColumnName = \'Y\', \'x-section\', \'stat\'')
		for s in range(nSyst):
			if (Hsyst[m][0][s] != 0):
				f.write(',\'error%05i\'' % s)
		f.write('\n')
		f.write('   Reaction = \'Bla\'\n')
		f.write('   Percent = ')
		for s in range(Hsyst[m, 0, :].sum() + 1):
			if (doPercent):
				f.write('true, ')
			else:
				f.write('false, ')
		f.write('\n')
		f.write('&END\n')

		# Loop over data point
		for d in range(nData):
			sf = Mdata[m][d]
			if (doPercent):
				sf = 100
			if (Mdata[m][d] != 0):
				f.write('%5.1f ' % d)
				f.write('%8.3f ' % Mdata[m][d])
				f.write('%8.3f ' % (Mstat[m][d] * sf))
				# Loop over systematics
				for s in range(nSyst):
					if (Hsyst[m][0][s] != 0):
						f.write('%8.3f ' % (Msyst[m][d][s] * sf))
				f.write('\n')
		f.write('\n')
		f.close()


# Dataset generator is a script to produce toy datasets in order to test averager
def genData(nData, nMes, nSyst, zMes, zSyst, doPoisson, doSmear, doPercent, seed, Stat, Syst):

	# Seed:
	np.random.seed(seed)

	# Parameters of uncertainties
	# Minumal and maximal relative statistical/systematic uncertainty
	vStat = np.array((Stat).split('-'), dtype=float)*0.01
	vSyst = np.array((Syst).split('-'), dtype=float)*0.01

	print 'Range of statistical uncertainties:', vStat
	print 'Range of systematic uncertainties:', vSyst

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

	return Mdata, Msyst, Mstat, Hsyst, Tshift, Tdata


def run():
	########################
	# setup input parameters
	########################
	usage = "DatasetGen.py parameters \nCreate toy dataset with gaussian uncertainties.\n"
	parser = OptionParser(usage)
	parser.add_option("-d", dest="nData", default=50, type="int", help="Number of data points")
	parser.add_option("-m", dest="nMes", default=2, type="int", help="Number of measurements")
	parser.add_option("-s", dest="nSyst", default=10, type="int", help="Number of systematic sources")

	parser.add_option("-M", dest="zMes", default=0.99, type="float", help="Fraction of non-emply measurements")
	parser.add_option("-S", dest="zSyst", default=0.99, type="float", help="Fraction of non-emply systematics")

	parser.add_option("--Percent", "-p", action="store_true", dest="percent", default=False,
					  help="True: relative uncertainty in %, False: Absolute uncertainty")
	parser.add_option("--Poisson", action="store_true", dest="poisson", default=False,
					  help="Use poissonian statistical uncertainties")
	parser.add_option("--Smear", action="store_true", dest="smear", default=False,
					  help="Systematic uncertainties are affected by statistics")
	parser.add_option("--Seed", dest="Seed", default=None, type="int",
					  help="Specify seed (default: take system time)")

	parser.add_option("--Stat", dest="Stat", default='2.5-4.5', type="string",
					  help="Range of stat uncertainty given in percent")
	parser.add_option("--Syst", dest="Syst", default='3.0-8.0', type="string",
					  help="Range of syst uncertainty given in percent")

	parser.add_option("--Fortran", "-f", action="store_true", dest="fortran", default=False,
					  help="True: Write output for fortran version of the averager")

	parser.add_option("--Plot", action="store_true", dest="plot", default=False,
					  help="True: Plot correlation matrix")

	parser.add_option("--Truth", "-t", action="store_true", dest="truth", default=True,
					  help="True: Store 'truth' information")

	for option in parser.option_list:
		if option.default != ("NO", "DEFAULT"):
			option.help += (" " if option.help else "") + "[default: %default]"

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
	doPercent = options.percent

	# random seed
	seed = options.Seed

	# size of uncertainties
	Stat = options.Stat
	Syst = options.Syst

	Mdata, Msyst, Mstat, Hsyst, Tshift, Tdata = genData(nData, nMes, nSyst, zMes, zSyst,
				doPoisson, doSmear, doPercent, seed, Stat, Syst)

	# Sture truth information
	if options.truth:
		np.savetxt('Tshift.out', Tshift, fmt='%1.3f')
		np.savetxt('Tdata.out', Tdata, fmt='%1.3f')

	# Write output files. Python format
	WriteCSV(Mdata, Msyst*Hsyst, Mstat)

	# Write output files. Fortran format
	if options.fortran:
		writeDAT(Mdata, Msyst, Mstat, Hsyst, doPercent)

	# Plot correlation matrix
	if options.plot:
		Corr = AveUtils.GetCorrMatrix(Msyst.reshape(nData*nMes,nSyst))

		AvePlot.PlotMatrix(Corr, fname='CorrGen.pdf',
						   xtitle='bin number', ytitle='bin number', label='Correlation')

if __name__ == "__main__":
	run()
