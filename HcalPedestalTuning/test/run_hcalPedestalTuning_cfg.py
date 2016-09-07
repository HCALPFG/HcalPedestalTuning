import FWCore.ParameterSet.Config as cms

process = cms.Process("USER")

process.load('Configuration.StandardSequences.Services_cff')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.default.limit = 10
process.MessageLogger.cerr.FwkReport.reportEvery = 1 

process.load('Configuration.EventContent.EventContent_cff') 
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.Geometry.GeometryDB_cff')
process.load("CondCore.DBCommon.CondDBSetup_cfi")

process.hcal_db_producer = cms.ESProducer("HcalDbProducer",
   dump = cms.untracked.vstring(''),
   file = cms.untracked.string('')
)

# Global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_v10', '')

# Specifiy an emap if you want to use something other than in GT
#process.es_pool = cms.ESSource("PoolDBESSource",
#        process.CondDBSetup,
#        timetype = cms.string('runnumber'),
#        toGet = cms.VPSet(
#            cms.PSet(
#                record = cms.string(
#                    "HcalElectronicsMapRcd"
#                    ),
#                tag = cms.string(
#                    "HcalElectronicsMap_v8.0_hlt"
#                    )
#                )
#            ),
#        connect = cms.string(
#            'frontier://FrontierProd/CMS_CONDITIONS'),
#            authenticationMethod = cms.untracked.uint32(0)
#        )    
#process.es_prefer_es_pool = cms.ESPrefer('PoolDBESSource', 'es_pool')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("HcalTBSource",
    fileNames = cms.untracked.vstring(
        # v1: input pedestals are default
        #"root://eoscms//eos/cms/store/group/dpg_hcal/comm_hcal/LS1/USC_272944.root"
        "file:/afs/cern.ch/user/j/jaehyeok/public/USC_274451.root"
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
    #ExceptionEmptyData = cms.untracked.bool(False),
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
    mapIOV             = cms.int32(7),
    tunePedestals      = cms.bool(True),
    addChannels        = cms.bool(True),
    pathToOldXMLBricks = cms.untracked.string('OldXML'),
    pathToNewXMLBricks = cms.untracked.string('NewXML'),
    pathToFinalXMLBricks = cms.untracked.string('FinalXML'),
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
    tagName            = cms.string('2016-may-18_HCAL_Pedestals')
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('output_histo.root')
)

#Paths
process.raw2digi_step = cms.Path(process.hcalDigis)
process.tuning_step = cms.Path(process.pedTuner)

# Schedule
process.schedule = cms.Schedule(process.raw2digi_step,process.tuning_step)

