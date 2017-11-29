from numpy import *
from numpy import linalg as LA

# Reduce correlation model
# shape of array syst: 
# [number of systematics, number of bins, number of measurements]
def ReduceCorrelationModel(syst, stat=[0], coef=0):

	nSyst = syst.shape[0]
	nBins = syst.shape[1]
	nMes = syst.shape[2]

	# Reshape array of systematics to
	# [number of bins * number of measurements, number of systematics]
	syst = syst.transpose().reshape(nBins*nMes,nSyst)
	stat = stat.transpose().ravel()

	# Calculate covariance matrix
	Cov = syst.dot(swapaxes(syst,0,1))
	if(len(stat)>1):
		Cov = Cov + diag(stat**2)

	print 'Cov \n', Cov

	# Decompose it in terms of nuisance parameters
	w, v = LA.eigh(Cov)
	# Avoid too small values	
	w = w*(w>1e-14)
	wNorm = sqrt(w)/(sum(sqrt(w)))
	print 'Eig \n', w
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

	print 'New total error: \n', sqrt((v**2).sum(axis=0))
	

	# convert to original shape of systematics array
	syst = NewSyst.reshape(nMes,nBins,nSystNew).transpose()
	stat = NewStat.reshape(nMes,nBins).transpose()
	return syst, stat

