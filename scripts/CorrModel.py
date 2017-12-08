from numpy import *
from numpy import linalg as LA

# Reduce correlation model
# shape of array syst: 
# [number of systematics, number of bins, number of measurements]
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

