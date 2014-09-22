
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
