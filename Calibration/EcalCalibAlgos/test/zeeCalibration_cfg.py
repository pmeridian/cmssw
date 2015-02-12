import FWCore.ParameterSet.Config as cms

###### ------------   configuration --------------------- #######
isMC = True
#prescale = 1

# electron cuts  
useGolden = False     
ELECTRON_ET_CUT_MIN = 20.0
ELECTRON_CUTS = "(abs(superCluster.eta)<2.5) && (ecalEnergy*sin(superClusterPosition.theta)>" + str(ELECTRON_ET_CUT_MIN) + ")"

# to read the HLT menu - should be always this one
HLTProcessName = "HLT"

# mass cuts (for T&P)
MASS_CUT_MIN = 60.

###### ------------   configuration --------------------- #######


processName='ZEECALIB'
process = cms.Process(processName)

#process.prescaler = cms.EDFilter("Prescaler",
#                                    prescaleFactor = cms.int32(prescale),
#                                    prescaleOffset = cms.int32(0)
#                                    )


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
#readFiles = cms.untracked.vstring()
#readFiles.extend( [
#        "file:zeeSkimChiara.root"
#        ] )
#process.source = cms.Source("PoolSource",
#                            fileNames = readFiles,
#)

if (not isMC):
    process.source.lumisToProcess = goodLumis                            


# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
if (not isMC):
    process.GlobalTag.globaltag = 'GR_R_42_V17::All'    #chiara, ancora da cambiare
else:
    process.GlobalTag = GlobalTag(process.GlobalTag, 'PHYS14_25_V1', '') 


# chiara: dovrebbe funzionare ancora con la nuova reco, ma controlla
# chiara: dove prende il file con le calibrazioni?
process.load('RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi')
process.load("RecoLocalCalo.EcalRecProducers.ecalRecalibRecHit_cfi")
process.ecalRecHit.doEnergyScale = cms.bool(False)       
process.ecalRecHit.doIntercalib = cms.bool(True)
process.ecalRecHit.doLaserCorrection = cms.bool(False)
process.ecalRecHit.EBRecHitCollection = "reducedEcalRecHitsEB"          
process.ecalRecHit.EERecHitCollection = "reducedEcalRecHitsEE"          
process.ecalRecHit.EBRecalibRecHitCollection = "EcalRecHitsEB"          # chiara: si chiamano EcalRecHitsEB ma sono reduced. Da capire 
process.ecalRecHit.EERecalibRecHitCollection = "EcalRecHitsEE"          # quale formato ci sara'. Inoltre, considera che il PF dopo
                                                                        # avrebbe bisogno i non reduced! che sono usati in
                                                                        # RecoParticleFlow/PFClusterProducer/interface/PFEcalRecHitCreator.h

# chiara:
# non so come cambiare il nome alla collezione di rechits reduced che deve diventare EcalRecHitsES invece di reducedEcalRecHitsES         
# giro questo modulo dummy che si limita a cambiare il nome.                                                                      
# non servira' quando girero' sui raw => NB: a me servono i rechits tutti, ma negli AOD ci sono solo i reducedEcalRecHitsES            
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

# chiara: bisogna spegnere la regression? 
process.particleFlowSuperClusterECAL.PFBasicClusterCollectionBarrel = cms.string('recalibParticleFlowBasicClusterECALBarrel')
process.particleFlowSuperClusterECAL.PFSuperClusterCollectionBarrel = cms.string('recalibParticleFlowSuperClusterECALBarrel')
process.particleFlowSuperClusterECAL.PFBasicClusterCollectionEndcap = cms.string('recalibParticleFlowBasicClusterECALEndcap')
process.particleFlowSuperClusterECAL.PFSuperClusterCollectionEndcap = cms.string('recalibParticleFlowSuperClusterECALEndcap')
process.particleFlowSuperClusterECAL.PFBasicClusterCollectionPreshower = cms.string('recalibParticleFlowBasicClusterECALPreshower')
process.particleFlowSuperClusterECAL.PFSuperClusterCollectionEndcapWithPreshower = cms.string('recalibParticleFlowSuperClusterECALEndcapWithPreshower')

# chiara: NB: in particleFlowSuperClusterECAL_cfi cambio solo i nomi di output e non di input. 
# chiara: NB: non ho cambiato i nomi delle associazioni in output

## new SCs - electrons association
## salva nell'evento gli elettroni e i SC associati ma solo in EE. Perche'?
#process.load("Calibration.EcalCalibAlgos.gedElectronRecalibSCAssociator_cfi")
#process.gedElectronRecalibSCAssociator.ebSuperclusters = cms.string('recalibParticleFlowSuperClusterECALBarrel')     
#process.gedElectronRecalibSCAssociator.eeSuperclusters = cms.string('recalibParticleFlowSuperClusterECALEndcap')     
#process.gedElectronRecalibSCAssociator.gsfElectrons = 'PassingWP90'        


# chiara: algo di calibrazione vero e proprio
process.load("Calibration.EcalCalibAlgos.zeeCalibration_cff")

# chiara: parametri eventualmente da cambiare
process.looper.maxLoops = cms.untracked.uint32(1)              # chiara: era 7
process.looper.electronSelection = cms.untracked.int32(-1)     # 0-1-2-3-4; -1 to do nothing
process.looper.histoFile = cms.string('myHistograms_test.root')
process.looper.zeeFile = cms.string('myZeePlots_test.root')
process.looper.initialMiscalibrationBarrel = cms.untracked.string('')
process.looper.initialMiscalibrationEndcap = cms.untracked.string('')
process.looper.ZCalib_CalibType = cms.untracked.string('RING')
process.looper.ZCalib_InvMass = cms.untracked.string('SCTRMass')
# chiara: dovrebbero restare cosi'
process.looper.rechitProducer   = cms.string('ecalRecHit')                          
process.looper.rechitCollection = cms.string('EcalRecHitsEB') 
process.looper.erechitProducer  = cms.string('ecalRecHit')                          
process.looper.erechitCollection = cms.string('EcalRecHitsEE')
process.looper.ebSuperclusters = cms.InputTag("particleFlowSuperClusterECAL","recalibParticleFlowSuperClusterECALBarrel",processName)
process.looper.eeSuperclusters = cms.InputTag("particleFlowSuperClusterECAL","recalibParticleFlowSuperClusterECALEndcapWithPreshower",processName)
process.looper.electrons = cms.InputTag("gedGsfElectrons","","RECO")
#process.looper.electronProducer = cms.string('PassingWP90')    # chiara: se spengo il filtro questo va spento
process.looper.HLTriggerResults = cms.InputTag("TriggerResults","",HLTProcessName)    
#Setting to null value avoids reading mc infos
#process.looper.mcProducer = cms.untracked.string('')           

## -----------------------------------------------------------
# electronID 
process.goodElectrons = cms.EDFilter("GsfElectronRefSelector",     
                                 src = cms.InputTag( 'gedGsfElectrons' ),
                                 cut = cms.string( ELECTRON_CUTS )             
                             )

process.PassingWP90 = process.goodElectrons.clone(
    cut = cms.string(
        process.goodElectrons.cut.value() +
        " && (ecalEnergy*sin(superClusterPosition.theta)>" + str(ELECTRON_ET_CUT_MIN) + ")"
        " && gsfTrack.isAvailable()"
        " && gsfTrack.hitPattern().numberOfLostHits(\'MISSING_INNER_HITS\') < 2"
        " && ( pfIsolationVariables().sumChargedHadronPt + pfIsolationVariables().sumNeutralHadronEt + pfIsolationVariables().sumPhotonEt )/pt < 0.25"
        " && abs(1./ecalEnergy - 1./(gsfTrack.p()))<0.05"   
        " && ( (isEB"
        " && (sigmaIetaIeta<0.01)"
        " && ( abs(deltaPhiSuperClusterTrackAtVtx)<0.15 )"
        " && ( abs(deltaEtaSuperClusterTrackAtVtx)<0.007 )"
        " && (hadronicOverEm<0.12)"
        ")"
        " || ( isEE"
        " && (sigmaIetaIeta<0.03)"
        " && ( abs(deltaPhiSuperClusterTrackAtVtx)<0.10 )"
        " && ( abs(deltaEtaSuperClusterTrackAtVtx)<0.009 )"
        " && (hadronicOverEm<0.1) "
        ") )"
        )
    )

# chiara, questo e' ancora vero?
# oltre a questo c'e' il parametro "electronSelection" che seleziona sulla classe
if (useGolden):
    process.PassingWP90.cut = process.PassingWP90.cut.value() + "&& (classification != 0)"   

process.Zele_sequence = cms.Sequence(
    process.goodElectrons * process.PassingWP90
    )


# chiara, ancora da rivedere
import copy
from HLTrigger.HLTfilters.hltHighLevel_cfi import *
process.ZEEHltFilter = copy.deepcopy(hltHighLevel)
process.ZEEHltFilter.throw = cms.bool(False)
process.ZEEHltFilter.HLTPaths = ["HLT_Ele*"]

# pairs
process.tagGsf =  cms.EDProducer("CandViewShallowCloneCombiner",
                                 decay = cms.string("PassingWP90 PassingWP90"),        
                                 checkCharge = cms.bool(False),
                                 cut   = cms.string("mass > " + str(MASS_CUT_MIN))
            )

process.tagGsfCounter = cms.EDFilter("CandViewCountFilter",
                                     src = cms.InputTag("tagGsf"),
                                     minNumber = cms.uint32(1)
                                     )

process.tagGsfFilter = cms.Sequence(process.tagGsf * process.tagGsfCounter)

process.tagGsfSeq = cms.Sequence(  )
if (not isMC):
    process.tagGsfSeq *= ( process.ZEEHltFilter )

process.tagGsfSeq *=  (  process.Zele_sequence * process.tagGsfFilter )

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string("/tmp/crovelli/testOut.root")
                               )

#process.Timing = cms.Service("Timing")

process.zFilterPath = cms.Path( # process.tagGsfSeq *     #chiara: capire perche' se lo accendo e non trova niente continua lo stesso e crasha dopo
                                                          # per ora lo spengo, ma va capito
                                process.ecalRecHit * process.ecalPreshowerRecHit * 
                                process.particleFlowRecHitPS * process.particleFlowClusterPS *    
                                process.particleFlowRecHitECAL * process.particleFlowClusterECALUncorrected * process.particleFlowClusterECAL *
                                process.particleFlowSuperClusterECAL )
                                ## * process.gedElectronRecalibSCAssociator )

#process.outpath = cms.EndPath(process.out)
