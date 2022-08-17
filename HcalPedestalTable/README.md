# HcalPedestalTuning/HcalPedestalTable 

`HCALPedestalTableMaker.C` is the code that makes pedestal table off PFGntuple. The usage is  
```
$ root -b -q HCALPedestalTableMaker.C++  
```
You should make a PFGntuple using a local pedestal run first. The following are two useful scripts for analysis of the pedestal measurement.

Script that compares two pedestal tables and generate a root file with histrogams of comparison:
```
$ python ped_compare.py [table1.txt] [table2.txt] [output.root]
``` 

Script that visualized a pedestal table or the result of comparison (output file of `ped_compare.py`):
```
$ root -b -q HCALPedestalAnalysis.C++  
```

After generating a pedestal table, check if there are missing channels. <br/> You can check it by scanning 2D plots, the output file of `HCALPedestalAnalysis.C++`. If there are too many missing channels to interpolate by hand, you can make a list of missing channels and run this code.
```
$ python interpolation.py [table.txt] [missing_channels.txt] [output.txt]
```
It returns the average means/width of surrounding channels as means/widths for these channels. The `missing_channels.txt` should follow the format below. 
```
#U ADC  << this is the unit
#             eta             phi             dep             det     cap0     cap1     cap2     cap3 widthcap0 widthcap1 widthcap2 widthcap3      DetId
                9              51               1              HB
               10              51               1              HB
                .               .               .               .
                .               .               .               .
                .               .               .               .                                   
```


Some notes:  

* Check if there are channels with large mean or rms. These need further investigation.  
* The default method used in `HCALPedestalTableMaker.C` is option 2. This uses the ADC bins from 0 to max+6 where max is the index of the peak. Therefore, the mean and the rms of the ADC histogram and the actual measurement can be sligtly different. 


