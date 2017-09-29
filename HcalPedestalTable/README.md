# HcalPedestalTuning/HcalPedestalTable 

`HCALPedestalTableMaker.C` is the code that makes pedestal table off PFGntuple. The usage is  
```
$ root -b -q HCALPedestalTableMaker.C++  
```
You should make a PFGntuple using a local pedestal run first. The following are two useful scripts for analysis of the pedestal measurement. 

Script that compares two pedestal tables and generate a root file with histrogams of comparison:
```
$ python ped_compare.py [table1] [table2] 
``` 

Script that visualized a pedestal table or the result of comparison (output file of `ped_compare.py`):
```
$ root -b -q HCALPedestalAnalysis.C++  
```

Some notes:  

* Check if there are channels with large mean or rms. These need further investigation.  
* After generating a pedestal table, check if there's missing channels. The means and the widths for these channels should be taken from the previous measurement, i.e., what is currently in the DB.
* The default method used in `HCALPedestalTableMaker.C` is option 2. This uses the ADC bins from 0 to max+6 where max is the index of the peak. Therefore, the mean and the rms of the ADC histogram and the actual measurement can be sligtly different. 


