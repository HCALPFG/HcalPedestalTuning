
Instructions on how to do run the pedestal tuning code 
=== 

**README is still in progress**

### Overview

The procedure of pedestal tuning is 
```
(1) take a pedestal run
(2) check the pedestal level of all channels. If they are close to the target value, we are done. If not, go to (3). 
(3) run this package to generate new DAC setting (xml files)
(4) copy the new setting to somewhere on cmsusr machine 
(5 modify the pedestal run snippet to read the new setting
(6) go to (1)
```
These steps are repeated untill all the channels have the desired pedestal level. 

This package takes care of the step (3); it extracts the the pedestal 
info from a local pedestal run, compare it with the target pedestal level, 
adjust the DAC setting, and produce the new XML files.    

### Structure of the code 

* readSettings
    * Read the current setting from xml files in Old_HBHEHF (let's change this name)

* tunePedestals
    * Loop over channelList , compare the pedestals, and set the new setting
    * The list of tuned channel is tunedChannelsByRBX

* writeOutput
    * Writes XML bricks
        * Makes a list of channels by RBX name
    * For each RBX, write a BRICK
    * Writes a txt files with all channels
    * writeUnstable
        * Writes a text file for the channels with
    * writeTuned
        * Writes a text file for the channels that are tuned based on tunedChannelsByRBX
    * addChannels 
        * Add channels that do not exist in the new emap, but do in the old emap

## What to change  

There are only a few places you need to change:

## How to run 

```
$cmsRun test/run_hcalPedestalTuning_cfg.py
```


## Reference  

[1] https://twiki.cern.ch/twiki/bin/viewauth/CMS/PedestalTuning

