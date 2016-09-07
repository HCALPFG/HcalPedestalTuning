
Instructions on how to do run the pedestal tuning code 
=== 

### Overview

The procedure of pedestal tuning is 
```
(1) take a local pedestal run 
(2) run this package to generate new LV setting
(3) put the new setting in the DB
```
These steps are repeated untill all the channels have the desired pedestal level. 

This package takes care of the step (2); it extracts the the pedestal 
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


Create a tag; this is necessary for the bricks, once they are uploaded
to the settings database, to be found by someone configuring a run.
Create a set of bricks with all settings 4.  Go into the code
NewPedestalTuning.cc and on Line 588, (temporarily) change DAC to “4”.
In the same code, on line 527, change tagName_ to “whatever name you
want” To change the label, go to newpedestal_cfg.py and on line 71,
put in a new label for “test_bricks_new”
 
You then have to send the “4” bricks to the settings database.  The
instructions for doing that are at:

https://twiki.cern.ch/twiki/bin/view/CMS/WebHome?topic=OnlineHCALDataSub
missionProceduresTOProdOMDSP5Server 

But the idea is this: zip up all of your xml files into a .zip file
(zip –r zipfile.zip directoryName).
 
Get the CaloOnlineTools/HcalOnlineDb package from CVS: 
cmsrel <CMSSW_X_X_X> (this gets the CMSSW package if you don’t already have it) 
Cd <CMSSW_X_X_X/src> 
cmsenv 
addpkg CaloOnlineTools/HcalOnlineDb 
scram b 
cd CaloOnlineTools/HcalOnlineDc/test 
gmake 
 
For help, type: ./xmlToolsRun 
To run it (have the name of the directory with your .xml bricks ready), 
type: ./configDbXml.sh ad follow the instructions. 
 
Load them into the database: 
 
/home/apeterma/PEDS/CMSSW_5_0_0/src/NewPedTuner/NewPedestalTuning/Fe
b26_2012_AllHcalTunedPeds 
