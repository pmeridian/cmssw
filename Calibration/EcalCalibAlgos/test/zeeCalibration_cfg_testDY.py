import FWCore.ParameterSet.Config as cms

###### ------------   configuration --------------------- #######
isMC = True
#prescale = 1

# to read the HLT menu - should be always this one
HLTProcessName = "HLT"

###### ------------   configuration --------------------- #######

processName='ZEECALIB'
process = cms.Process(processName)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')  
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# chiara
from Calibration.EcalCalibAlgos.DYJetsToLLPhys14SkimZee_cff import *
process.source = cms.Source("PoolSource",
                            fileNames = readFiles,
)

if (not isMC):
    process.source.lumisToProcess = goodLumis                            

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
if (not isMC):
    process.GlobalTag.globaltag = 'GR_R_42_V17::All'    #chiara, ancora da cambiare
else:
    process.GlobalTag = GlobalTag(process.GlobalTag, 'PHYS14_25_V1', '') 


process.load('RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi')
process.load("RecoLocalCalo.EcalRecProducers.ecalRecalibRecHit_cfi")
process.ecalRecHit.doEnergyScale = cms.bool(False)       
process.ecalRecHit.doIntercalib = cms.bool(True)
process.ecalRecHit.doLaserCorrection = cms.bool(False)
process.ecalRecHit.EBRecHitCollection = "reducedEcalRecHitsEB"          
process.ecalRecHit.EERecHitCollection = "reducedEcalRecHitsEE"          
process.ecalRecHit.EBRecalibRecHitCollection = "EcalRecHitsEB"          
process.ecalRecHit.EERecalibRecHitCollection = "EcalRecHitsEE"          
                                                                        
process.load("RecoLocalCalo.EcalRecProducers.esDummyRecHit_cfi")
process.ecalPreshowerRecHit.ESRecHitCollection = "reducedEcalRecHitsES"

# PF rechits
process.load('RecoParticleFlow.PFClusterProducer.particleFlowRecHitECAL_cfi')
process.load('RecoParticleFlow.PFClusterProducer.particleFlowRecHitPS_cfi') 

# PF clustering       
process.load('RecoParticleFlow.PFClusterProducer.particleFlowClusterECALUncorrected_cfi')
process.load('RecoParticleFlow.PFClusterProducer.particleFlowClusterECAL_cfi')
process.load('RecoParticleFlow.PFClusterProducer.particleFlowClusterPS_cfi') 
process.load("RecoEcal.EgammaClusterProducers.particleFlowSuperClusterECAL_cfi")

process.particleFlowSuperClusterECAL.PFBasicClusterCollectionBarrel = cms.string('recalibParticleFlowBasicClusterECALBarrel')
process.particleFlowSuperClusterECAL.PFSuperClusterCollectionBarrel = cms.string('recalibParticleFlowSuperClusterECALBarrel')
process.particleFlowSuperClusterECAL.PFBasicClusterCollectionEndcap = cms.string('recalibParticleFlowBasicClusterECALEndcap')
process.particleFlowSuperClusterECAL.PFSuperClusterCollectionEndcap = cms.string('recalibParticleFlowSuperClusterECALEndcap')
process.particleFlowSuperClusterECAL.PFBasicClusterCollectionPreshower = cms.string('recalibParticleFlowBasicClusterECALPreshower')
process.particleFlowSuperClusterECAL.PFSuperClusterCollectionEndcapWithPreshower = cms.string('recalibParticleFlowSuperClusterECALEndcapWithPreshower')

process.load("Calibration.EcalCalibAlgos.zeeCalibration_cff")
process.looper.maxLoops = cms.untracked.uint32(3)              # chiara: era 7
process.looper.electronSelection = cms.untracked.int32(-1)     # 0-1-2-3-4; -1 to do nothing
process.looper.histoFile = cms.string('myHistograms_test.root')
process.looper.zeeFile = cms.string('myZeePlots_test.root')
process.looper.initialMiscalibrationBarrel = cms.untracked.string('')
process.looper.initialMiscalibrationEndcap = cms.untracked.string('')
process.looper.ZCalib_CalibType = cms.untracked.string('RING')
process.looper.ZCalib_InvMass = cms.untracked.string('SCTRMass')
process.looper.rechitProducer   = cms.string('ecalRecHit')                          
process.looper.rechitCollection = cms.string('EcalRecHitsEB') 
process.looper.erechitProducer  = cms.string('ecalRecHit')                          
process.looper.erechitCollection = cms.string('EcalRecHitsEE')
process.looper.ebSuperclusters = cms.InputTag("particleFlowSuperClusterECAL","recalibParticleFlowSuperClusterECALBarrel",processName)
process.looper.eeSuperclusters = cms.InputTag("particleFlowSuperClusterECAL","recalibParticleFlowSuperClusterECALEndcapWithPreshower",processName)
process.looper.electrons = cms.InputTag("gedGsfElectrons","","RECO")
process.looper.HLTriggerResults = cms.InputTag("TriggerResults","",HLTProcessName)    
#Setting to null value avoids reading mc infos
#process.looper.mcProducer = cms.untracked.string('')           

process.out = cms.OutputModule("PoolOutputModule",
                               outputCommands = "drop *_*_*_*",
                               fileName = cms.untracked.string("/tmp/crovelli/testOut.root")
                               )

process.zFilterPath = cms.Path( process.ecalRecHit * process.ecalPreshowerRecHit * 
                                process.particleFlowRecHitPS * process.particleFlowClusterPS *    
                                process.particleFlowRecHitECAL * process.particleFlowClusterECALUncorrected * process.particleFlowClusterECAL *
                                process.particleFlowSuperClusterECAL )

process.outpath = cms.EndPath(process.out)
