**Blas:**

Performance test of Haverager is demonstrated with python script, which generates random data samples and averaged them. HAverager with and without blas library is compared. 
```
./test.py
```
**Blue:**

This test performs a comparison of different ways of systematic evaluation of HAverager between each other and with Blue averager.

Comparion is performed using Zpt2012 ATLAS data. Blue averager is precomputed. Output of Blue averager is stored in files ZPtBlueCor.txt and ZptBlueOut.txt. Files myBLUEZptTest.cxx and myBLUEZptTest.inp have been used to perform a calculation with Blue.

Evaluation of combination with HAverager requires csv files of Zpt data (can be found here ../Zpt/*csv).
```
cp ../Zpt/*csv .
./testBlue.py
```
**Closure:**

This test performs a comparison of averaged data with "truth", taken from random number generator. In order to run this test random data sample have to be generated using script DatasetGen.py.
```
DatasetGen.py -d 10 -s 5
./test.py
```
**Zpt:**

This test performs a comparison of fortran and python implementation of the averager. Zpt2012 ATLAS data are usind. In order to run fortran version of the averager existing csv files have to be converted to dat files using supplied script.
```
ConvertCSVtoDat.py *.csv
averager steeringZpt
./test.py
```
**ZY:**

This script demonstrate fortran and python implementation of HAverager using ZY2011 ATLAS data.
```
ConvertDattoCSV.py [el,mu]*dat
averager steering
./test.py
```
