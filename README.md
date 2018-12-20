**Commands for installation:**  
```
autoreconf --install  
./configure  --enable-python  
make  
make install
```

**Usage:**  
```
source iniAverager  
cd ./test/Zpt  
```
**Create .dat files fron .csv (requires python pandas):**
```
ConvertCSVtoDat.py *.csv
averager steeringZpt  
plotAve.py TOutF
```
**Results**  
The results are printed out in the screen and stored (in case of default steering) in the directory ./TOutF
and illustrated in several plots in current directory

**More example can be found in the test folder**