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
process.GlobalTag.globaltag = "GR_P_V42_AN4::All"
process.es_prefer_GlobalTag = cms.ESPrefer('PoolDBESSource','GlobalTag')
process.prefer("GlobalTag")

### Read the latest HCAL Emap
#process.GlobalTag.toGet = cms.VPSet(
#    cms.PSet(record = cms.string("HcalElectronicsMapRcd"),
#             tag = cms.string("HcalElectronicsMap_v7.04_hlt"),
#             #connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_GLOBALTAG")
#             connect = cms.untracked.string("frontier://PromptProd/CMS_COND_31X_HCAL")
#             )
#    )
###

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("HcalTBSource",
    fileNames = cms.untracked.vstring(
#        'file:USC_213237.root'
#        'file:/afs/cern.ch/user/k/kkaadze/dataHOUSC/USC_216300.root'

      #Taken on 15 Jan                                
         'file:/afs/cern.ch/user/k/kkaadze/dataHOUSC/USC_217384.root',
         'file:/afs/cern.ch/user/k/kkaadze/dataHOUSC/USC_217384.1.root'
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
    debug              = cms.int32(1),  #needed for JS version of the code
    digisName          = cms.string("hcalDigis"),
    mapIOV             = cms.int32(5),
    tunePedestals      = cms.bool(True),
#    pathToOldXMLBricks = cms.untracked.string('Aug01_2013'),
    pathToOldXMLBricks = cms.untracked.string('Old_HOOnly'),
#    xmlBricksOutputDir = cms.untracked.string('test_bricks_new'),
    xmlBricksOutputDir = cms.untracked.string('New_HOOnly'),
    uniformDACSettings = cms.bool(False),
    initialDACSetting  = cms.int32(2),
#    SiPMRBXes          = cms.vstring('HO1P10', 'HO2P12'),
    SiPMRBXes          = cms.vstring('HO001', 'HO002', 'HO003', 'HO004', 'HO005', 'HO006',
                                     'HO007', 'HO008', 'HO009', 'HO010', 'HO011', 'HO012',
                                     'HO1M02', 'HO1M04', 'HO1M06', 'HO1M08', 'HO1M10', 'HO1M12',
                                     'HO1P02', 'HO1P04', 'HO1P06', 'HO1P08', 'HO1P10', 'HO1P12',
                                     'HO2M02', 'HO2M04', 'HO2M06', 'HO2M08', 'HO2M10', 'HO2M12',
                                     'HO2P02', 'HO2P04', 'HO2P06', 'HO2P08', 'HO2P10', 'HO2P12'),
    targetPedestalHPD  = cms.double(3.0),
    targetPedestalSiPM = cms.double(11.0),
    dac2adc            = cms.double(0.6),
    tagName            = cms.string('Aug08_2013')
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('alberto_histo.root')
)

#Paths
process.raw2digi_step = cms.Path(process.hcalDigis)
process.tuning_step = cms.Path(process.pedTuner)

# Schedule
process.schedule = cms.Schedule(process.raw2digi_step,process.tuning_step)

