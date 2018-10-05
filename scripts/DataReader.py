#!/usr/bin/env python

import numpy as np
import pandas as pd
import os
import glob
                
def hash_types(typee,list):                     # Internal function used to hash different error types and bins.
	for i in range(len(list)):
		#print i,typee,list[i]
		#print i,type(typee),type(list[i])
		if(type(list[i])==float):
			if list[i] == typee[0]:
				return i
		if list[i] == typee:
			return i
	raise AttributeError('No such type/bin')

def paverage(_bins,_data,_error):

# This routine parses all .csv data in the current directory for the Fortran routine average. 
#	_bins:  names of bin columns - input as 'bin1,bin2,...'
#	_data:  name of data column - input as 'data1'
#	_error: names of error columns - every column can have a modifier, note that one must have the Stat modifier, for statistical error.
#           - input as 'error1:corr,error2,error3:Stat,...'. Note different columns should have different names.

	# Define variables visible outsige of modile
	# Name of input files
	global fnames
	# Name of systematic uncertainties
	global oerror
	# Binning
	global bins
	# Bin names
	global binnames

	# Get all .csv in current directory
	path = os.getcwd()
	files = glob.glob(path + "/*.csv")
	fnames = glob.glob("*.csv")
	#files = glob.glob(path + "/*.txt")
	#fnames = glob.glob("*.txt")


	# Read all .csv in current directory and make a list of dataframes
	listofdfs = []
	for fileloc in files:
		print fileloc
		file = pd.read_csv(fileloc)
		#print file
		listofdfs += [file]

	# Find all the bin names in all the files by making a super dataframe 
	# containing all the individual dataframes, and then using the "groupby" function
	df = pd.DataFrame()
	for tempdf in listofdfs:
		df = df.append(tempdf)

	# Extruct the list of fields
	fields = df.columns.values.tolist()

	# Check input parameters. If they are not given, use default
	if(_data==''):
		data = 'data'
	else:
		data = _data

	if(_bins==''):
		bins = [s for s in fields if "bin" in s]
	else:
		bins  = _bins.split(',')
 
	#print fields

	if(_error==''):
		serror = [s for s in fields if "stat" in s]
		oerror = [s for s in fields if ("stat" not in s) and ("bin" not in s) and ("data" not in s)]
	else:
		stmp = _error.split(',')
		serror = [s for s in stmp if "stat" in s]
		oerror = [s for s in stmp if "stat" not in s]	

	# Group the mega-dataframe by the bins and
	# extract the bin names from the grouped dataframes. 
	grouped = list(df.groupby(bins))                 

	binnames = [item[0] for item in grouped] 

	#print 'Bin names: ', binnames

	# If the bin names are tupples - i.e. if more than one column specifies the bin.
	if(type(binnames[0])==tuple):                        
		binnames = [list(row) for row in binnames]   # Convert tuples to lists.

	n = len(binnames)
	m = len(files)
	l = len(oerror)

# Loop over all the dataframes and iteratively create lists of data and errors.
	data_ = np.zeros((n,m))
	serror_ = np.zeros((n,m))   
	oerror_ = np.zeros((l,n,m))     
	
	i = 0
	# Loop over dataframes
	for cdf in listofdfs:                              
		row_iterator = cdf.iterrows()
		#print 'File, row',row_iterator,cdf
		for index,row in row_iterator:                      # Over every row of every dataframe in listofdfs
			#print 'Row',row,index
			temp = row[bins]
			#print 'bin info:',temp                                # Find the bin
			cbin = temp.values.tolist()
			#print 'Cbin:',type(cbin),type(cbin[0]),cbin[0]
			hbin = hash_types(cbin,binnames)                # Hash it
			
			# Input Data
			tdata = row[data]                              # Extract the data
			data_[hbin,i] = tdata                          # And insert it

			# Input 'Stat' errors
			#print 'Find stat ',serror,row[serror]  
			tserror = row[serror]                         
			tserror1 = np.array(tserror.values.tolist()) # Parse tserror into a numpy array
			tserror2 = np.nan_to_num(tserror1)           # Clean up any possible "NaN" induced by missing values for error columns (i.e. if error is not present)
			tserror3 = np.linalg.norm(tserror2)          # Add all the values in quadrature
			#print tserror1,tserror2,tserror3
			serror_[hbin,i] = tserror3
			
			# Input other errors
			j = 0
			for coerror in oerror:                      # Read coerror as "current other error"
				#print coerror
				if coerror in row:
					toerror = row[coerror]
					oerror_[j,hbin,i] = toerror
				j += 1
		i += 1

	return data_,serror_,oerror_


