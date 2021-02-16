from numpy import *
from numpy import linalg as LA

# Reduce correlation model
# Script construct covariance matrix and
# decompose it to eigenvectors and eigenvalues.
# Eigenvalues are ordered from smallest to biggest.
# New covariance matrix is constructed using large
# eigenvalues and corresponding eigenvectors, while small
# are ignored. Sum of all ignored normalized eigenvalues are
# less then given parameter coef.
# shape of array syst: 
# [number of systematics, number of bins, number of measurements]
# stat uncertainty can be optionally added to diagonal of the
# covariance matrix.
# coef define a maximum of sum of ignored normalized eigenvalues.
def ReduceCorrelationModel(syst, stat=[0], coef=0):

	print 'Input shape ', syst.shape

	nSyst = syst.shape[0]
	nBins = syst.shape[1]

	if(len(syst.shape)==3):
		nMes = syst.shape[2]
	else:
		nMes = 1

	# Reshape array of systematics to
	# [number of bins * number of measurements, number of systematics]
	syst = syst.transpose().reshape(nBins*nMes,nSyst)

	# Calculate covariance matrix
	Cov = syst.dot(swapaxes(syst,0,1))
	if(len(stat)>1):
		stat = stat.transpose().ravel()
		Cov = Cov + diag(stat**2)

	# Decompose it in terms of nuisance parameters
	w, v = LA.eigh(Cov)
	# Avoid too small values	
	w = w*(w>1e-14)
	wNorm = sqrt(w)/(sum(sqrt(w)))
	#print 'Eig \n', w
	print 'Eig norm \n', wNorm
	# Count only big eigenvalues
	tmp=0
	counter=0
	for x in wNorm:
		tmp += x
		if(tmp>coef):
			break
		counter += 1;

	nSystNew = nMes*nBins - counter
	print 'Number of new systematic sources: ',nSystNew

	# Create array of new systematics
	v = v.dot(diag(sqrt(w)))

	# Reduce correlation model
	NewSyst = v[:,counter:]
 	NewStat = sqrt((v[:,:counter]**2).sum(axis=1))
	
	if(nMes > 1):
		# convert to original shape of systematics array
		syst = NewSyst.reshape(nMes,nBins,nSystNew).transpose()
		stat = NewStat.reshape(nMes,nBins).transpose()
	else:
		syst = NewSyst.transpose()
		stat = NewStat.transpose()
	return syst, stat

# Get correlation matrix from array of systematics
def GetCorrMatrix(A):
    Cov = GetCov(A)
    Corr = GetCorr(Cov)
    return Corr

# Get correlation from covariance matrix
def GetCorr(Cov):
	B  = diag(Cov)
	nB = B.size
	C = ones((nB,nB))
	D = B.reshape(nB,1)*C*B
	return Cov/sqrt(D)

# Get covariance matrix from array of systematics
def GetCov(A, axis=0):
	if(axis==0):
		Cov = A.dot(A.transpose())
	else:
		Cov = (A.transpose()).dot(A)
	return Cov