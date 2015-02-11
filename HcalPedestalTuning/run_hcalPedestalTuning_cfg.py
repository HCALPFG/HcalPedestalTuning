import FWCore.ParameterSet.Config as cms

process = cms.Process("USER")

process.load('Configuration.StandardSequences.Services_cff')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.default.limit = 10
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.Geometry.GeometryDB_cff')

process.hcal_db_producer = cms.ESProducer("HcalDbProducer",
   dump = cms.untracked.vstring(''),
   file = cms.untracked.string('')
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "GR_P_V49::All"
process.es_prefer_GlobalTag = cms.ESPrefer('PoolDBESSource','GlobalTag')
process.prefer("GlobalTag")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )

process.source = cms.Source("HcalTBSource",
    fileNames = cms.untracked.vstring(
#This run was taken with DB pedestals
#       'file:/afs/cern.ch/user/k/kkaadze/dataHOUSC/USC_218117.root'

#This run was taken with v1 tuned pedestals
#       'file:/afs/cern.ch/user/k/kkaadze/dataHOUSC/USC_218287.root'

#This run was taken with v2 tuned pedestals
#       'file:/afs/cern.ch/user/k/kkaadze/dataHOUSC/USC_218363.root'

#This run was taken with v3 tuned pedestals
#       'file:/afs/cern.ch/user/k/kkaadze/dataHOUSC/USC_218465.root'

#This run was taken with v4 tuned pedestals
#       'file:/afs/cern.ch/user/k/kkaadze/dataHOUSC/USC_218469.root'

#This run was taken with v5 tuned pedestals
#       'file:/afs/cern.ch/user/k/kkaadze/dataHOUSC/USC_218503.root'

#This run was taken with v6 tuned pedestals
#       'file:/afs/cern.ch/user/k/kkaadze/dataHOUSC/USC_218518.root'                                                                
#28/03/2014
#This run was taken with v7 tuned pedestals -- HO/Peds
#       'file:/afs/cern.ch/user/k/kkaadze/dataHOUSC/USC_220222.root'

#This run was taken with v8 tuned pedestals -- HO/Peds
#       'file:/afs/cern.ch/user/k/kkaadze/dataHOUSC/USC_220232.root'

#This run was taken with v9 tuned pedestals -- HO/Peds
#       'file:/afs/cern.ch/user/k/kkaadze/dataHOUSC/USC_220244.root'                                                                                                
#This run was taken with v10 tuned pedestals -- HO/Peds
#       'file:/afs/cern.ch/user/k/kkaadze/dataHOUSC/USC_220262.root'                                                                                                
#This run was taken with v11 tuned pedestals -- HO/Peds
#       'file:/afs/cern.ch/user/k/kkaadze/dataHOUSC/USC_220263.root'                                                                                                
#This run was taken with v12 tuned pedestals -- HO/Peds
#       'file:/afs/cern.ch/user/k/kkaadze/dataHOUSC/USC_220418.root'

#This run was taken with v13 tuned pedestals -- HO/Peds
#       'file:/afs/cern.ch/user/k/kkaadze/dataHOUSC/USC_220425.root'

#This run was taken with v14 tuned pedestals -- HO/Peds
       'file:/afs/cern.ch/user/k/kkaadze/dataHOUSC/USC_220426.root'                                                                                                                                
   )
)

process.options = cms.untracked.PSet( 
    wantSummary = cms.untracked.bool(False) ) ## default is false

# HCAL digis
process.hcalDigis = cms.EDProducer("HcalRawToDigi",
    # Flag to enable unpacking of ZDC channels (default = false)
    UnpackZDC = cms.untracked.bool(False),
    # Optional filter to remove any digi with "data valid" off, "error" on, 
    # or capids not rotating
    FilterDataQuality = cms.bool(True),
    # Do not complain about missing FEDs
    ExceptionEmptyData = cms.untracked.bool(False),
    InputLabel = cms.InputTag("source"),
    # Use the defaults for FED numbers
    # Do not complain about missing FEDs
    ComplainEmptyData = cms.untracked.bool(False),
    # Flag to enable unpacking of calibration channels (default = false)
    UnpackCalib = cms.untracked.bool(True),
    # At most ten samples can be put into a digi, if there are more
    # than ten, firstSample and lastSample select which samples
    # will be copied to the digi
    firstSample = cms.int32(0),
    lastSample = cms.int32(9)
)

process.pedTuner = cms.EDAnalyzer('HcalPedestalTuning',
    debug              = cms.int32(0),  
    digisName          = cms.string("hcalDigis"),
    mapIOV             = cms.int32(5),
    tunePedestals      = cms.bool(True),
    pathToOldXMLBricks = cms.untracked.string('Old_HBHEHF'),
    xmlBricksOutputDir = cms.untracked.string('New_HBHEHF'),
    uniformDACSettings = cms.bool(False),
    initialDACSetting  = cms.int32(0),
    isHO = cms.bool(True),                                  
    SiPMRBXes          = cms.vstring('HO001', 'HO002', 'HO003', 'HO004', 'HO005', 'HO006',
                                     'HO007', 'HO008', 'HO009', 'HO010', 'HO011', 'HO012',
                                     'HO1M02', 'HO1M04', 'HO1M06', 'HO1M08', 'HO1M10', 'HO1M12',
                                     'HO1P02', 'HO1P04', 'HO1P06', 'HO1P08', 'HO1P10', 'HO1P12',
                                     'HO2M02', 'HO2M04', 'HO2M06', 'HO2M08', 'HO2M10', 'HO2M12',
                                     'HO2P02', 'HO2P04', 'HO2P06', 'HO2P08', 'HO2P10', 'HO2P12'),
    targetPedestalHPD  = cms.double(3.0),
    targetPedestalSiPM = cms.double(9.0),
    dac2adc            = cms.double(0.6),
    tagName            = cms.string('2014-feb-09_HCAL_Pedestals')
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('output_histo.root')
)

#Paths
process.raw2digi_step = cms.Path(process.hcalDigis)
process.tuning_step = cms.Path(process.pedTuner)

# Schedule
process.schedule = cms.Schedule(process.raw2digi_step,process.tuning_step)

