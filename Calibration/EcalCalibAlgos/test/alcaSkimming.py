import FWCore.ParameterSet.Config as cms

#
#    _____             __ _                        _   _
#   / ____|           / _(_)                      | | (_)
#   | |     ___  _ __ | |_ _  __ _ _   _ _ __ __ _| |_ _  ___  _ __
#   | |    / _ \| '_ \|  _| |/ _` | | | | '__/ _` | __| |/ _ \| '_ \
#   | |___| (_) | | | | | | | (_| | |_| | | | (_| | |_| | (_) | | | |
#    \_____\___/|_| |_|_| |_|\__, |\__,_|_|  \__,_|\__|_|\___/|_| |_|
#                             __/ |
#                            |___/

# From which kind of dataset are you starting from?
MC = True
runFromAOD  = True
runFromALCA = False

# Do you want to reReco ECAL RecHits (also from RAW if you have them) ? 
ECALRecalib = False
# if ECALRecalib, then set these ones below:
ECALFromRAW = False
ApplyInterCalib = False
ApplyLaser = False   #Switch to turn on/off application of LC. To be set to False when running from RAW and want to produce LC=1

# Do you want to filter events? 
HLTFilter = False
HLTPath = "HLT_Ele"
HLTProcessName = "HLT"

#electron cuts
ELECTRON_ET_CUT_MIN = 20.0
ELECTRON_CUTS = "(abs(superCluster.eta)<2.5) && (ecalEnergy*sin(superClusterPosition.theta)>" + str(ELECTRON_ET_CUT_MIN) + ")"

#mass cuts (for Zee selection)
MASS_CUT_MIN = 60.


#    _____  __             _             _         _
#   / ____|/ _|           | |           | |       | |
#   | |    | |_ __ _   ___| |_ __ _ _ __| |_ ___  | |__   ___ _ __ ___
#   | |    |  _/ _` | / __| __/ _` | '__| __/ __| | '_ \ / _ \ '__/ _ \
#   | |____| || (_| | \__ \ || (_| | |  | |_\__ \ | | | |  __/ | |  __/
#    \_____|_| \__, | |___/\__\__,_|_|   \__|___/ |_| |_|\___|_|  \___|
#               __/ |
#              |___/
   

if (not runFromALCA):
    processName = 'ALCASKIM'
else:
    processName ='ALCARERECO'

process = cms.Process(processName)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
#process.load('Configuration.EventContent.EventContent_cff')
#process.load('SimGeneral.MixingModule.mixNoPU_cfi') 
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')  
#process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff') 
#process.load('Configuration.StandardSequences.RawToDigi_cff')
#process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
#process.load('Configuration.StandardSequences.L1Reco_cff')
#process.load('Configuration.StandardSequences.Reconstruction_cff')
#process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.AlCaRecoStreams_cff')


process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(200)
)

#from Calibration.EcalCalibAlgos.DoubleElectron_Jul05_ALCAELECTRON_cff import *
#from Calibration.EcalCalibAlgos.Cert_160404_172802_cff import *

readFiles = cms.untracked.vstring()

readFiles.extend( [
        "file:fileAOD.root"   # MC DY file from phys14
] )

process.source = cms.Source("PoolSource",
                            fileNames = readFiles
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag                        
if (MC):
    process.GlobalTag = GlobalTag(process.GlobalTag, 'PHYS14_25_V1', '')     
else:
    process.GlobalTag.globaltag = 'GR_R_42_V17::All'   # chiara, da cambiare


##    _____ _           _                     ___    _
##   | ____| | ___  ___| |_ _ __ ___  _ __   |_ _|__| |
##   |  _| | |/ _ \/ __| __| '__/ _ \| '_ \   | |/ _` |
##   | |___| |  __/ (__| |_| | | (_) | | | |  | | (_| |
##   |_____|_|\___|\___|\__|_|  \___/|_| |_| |___\__,_|
##


process.selectedElectrons = cms.EDFilter("GsfElectronRefSelector",
                                 src = cms.InputTag( 'gedGsfElectrons' ),
                                 cut = cms.string( ELECTRON_CUTS )
                             )

process.PassingWP90 = process.selectedElectrons.clone(
    cut = cms.string(
        process.selectedElectrons.cut.value() +
        # chiara: qui metto 2012id, WP loose
        " && (ecalEnergy*sin(superClusterPosition.theta)>" + str(ELECTRON_ET_CUT_MIN) + ")"
        " && gsfTrack.isAvailable()" 
        " && gsfTrack.hitPattern().numberOfLostHits(\'MISSING_INNER_HITS\') < 2"    # chiara: in 72x mi pare essere diventato cosi'... 
        # chiara: match con conversioni da mettere offline
        # chiara: IP cuts da mettere offline
        " && ( pfIsolationVariables().sumChargedHadronPt + pfIsolationVariables().sumNeutralHadronEt + pfIsolationVariables().sumPhotonEt )/pt < 0.25" 
        # chiara: metto un taglio piu' loose e non applico correzioni per PU 
        " && abs(1./ecalEnergy - 1./p4.P)<0.05"
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

process.ele_sequence = cms.Sequence(
    process.PassingWP90
    )


###############################
# ECAL Recalibration
###############################
process.electronRecalib = cms.Sequence()

if (ECALRecalib):
    process.load("Calibration.EcalCalibAlgos.electronRecalibSCAssociatorSH_cfi")
            
    if (not ECALFromRAW):
        process.load("RecoLocalCalo.EcalRecProducers.ecalRecalibRecHit_cfi")
        process.ecalRecHit.doIntercalib = cms.bool(ApplyInterCalib)
        process.ecalRecHit.doLaserCorrection = cms.bool(ApplyLaser)
        if (runFromAOD):
            process.ecalRecHit.EBRecHitCollection = "reducedEcalRecHitsEB"
            process.ecalRecHit.EERecHitCollection = "reducedEcalRecHitsEE"
        elif (runFromALCA): 
            process.ecalRecHit.EBRecHitCollection = "alCaIsolatedElectrons:alcaBarrelHits"
            process.ecalRecHit.EERecHitCollection = "alCaIsolatedElectrons:alcaEndcapHits"

        process.ecalRecHit.EBRecalibRecHitCollection = "EcalRecHitsEB"
        process.ecalRecHit.EERecalibRecHitCollection = "EcalRecHitsEE"
        process.electronRecalib *= (process.ecalRecHit)
    else:
        # restarting from ECAL RAW to reconstruct amplitudes and energies
        process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
        process.load('RecoLocalCalo.Configuration.RecoLocalCalo_cff')
        process.ecalRecHit.laserCorrection=cms.bool(ApplyLaser)
        # no switch in standard recHit producer to apply new intercalibrations
        process.electronRecalib *= ( (process.ecalDigis+process.ecalPreshowerDigis) * process.ecalLocalRecoSequence)
        
    process.load("RecoEcal.Configuration.RecoEcal_cff")
    process.correctedHybridSuperClusters.corectedSuperClusterCollection = 'recalibSC'
    process.correctedMulti5x5SuperClustersWithPreshower.corectedSuperClusterCollection = 'endcapRecalibSC'

    if (runFromAOD):
        process.multi5x5SuperClustersWithPreshower.preshRecHitProducer = cms.InputTag("reducedEcalRecHitsES")
        process.multi5x5PreshowerClusterShape.preshRecHitProducer = cms.InputTag("reducedEcalRecHitsES")
    elif (runFromALCA):
        # chiara: controlla se li hanno rimessi, perche' al 26/11/14 erano scomparsi
        process.multi5x5SuperClustersWithPreshower.preshRecHitProducer = cms.InputTag("alCaIsolatedElectrons","alcaPreshowerHits")
        process.multi5x5PreshowerClusterShape.preshRecHitProducer = cms.InputTag("alCaIsolatedElectrons","alcaPreshowerHits")

    process.electronRecalibSCAssociatorSH.scIslandCollection = cms.string('endcapRecalibSC')
    process.electronRecalibSCAssociatorSH.scIslandProducer = cms.string('correctedMulti5x5SuperClustersWithPreshower')
    process.electronRecalibSCAssociatorSH.scProducer = cms.string('correctedHybridSuperClusters')
    process.electronRecalibSCAssociatorSH.scCollection = cms.string('recalibSC')
    process.electronRecalibSCAssociatorSH.electronProducer = 'gedGfElectrons'
    process.electronRecalib *= (process.hybridClusteringSequence* process.multi5x5ClusteringSequence * process.multi5x5PreshowerClusteringSequence * process.electronRecalibSCAssociatorSH)

if (runFromAOD):
    process.alCaIsolatedElectrons.esRecHitsLabel = cms.InputTag("reducedEcalRecHitsES")
elif (runFromALCA):
    # chiara: controlla se li hanno rimessi, perche' al 26/11/14 erano scomparsi  
    process.alCaIsolatedElectrons.esRecHitsLabel = cms.InputTag("alCaIsolatedElectrons","alcaPreshowerHits")

if (ECALRecalib):
    process.alCaIsolatedElectrons.ebRecHitsLabel = cms.InputTag("ecalRecHit:EcalRecHitsEB")
    process.alCaIsolatedElectrons.eeRecHitsLabel = cms.InputTag("ecalRecHit:EcalRecHitsEE")
    process.alCaIsolatedElectrons.electronLabel = cms.InputTag("electronRecalibSCAssociatorSH")        
elif (runFromAOD):
    process.alCaIsolatedElectrons.ebRecHitsLabel = cms.InputTag("reducedEcalRecHitsEB")
    process.alCaIsolatedElectrons.eeRecHitsLabel = cms.InputTag("reducedEcalRecHitsEE")
elif (runFromALCA):
    process.alCaIsolatedElectrons.ebRecHitsLabel = cms.InputTag("alCaIsolatedElectrons:alcaBarrelHits")
    process.alCaIsolatedElectrons.eeRecHitsLabel = cms.InputTag("alCaIsolatedElectrons:alcaEndcapHits")
    
##    ____       _
##   |  _ \ __ _(_)_ __ ___
##   | |_) / _` | | '__/ __|
##   |  __/ (_| | | |  \__ \
##   |_|   \__,_|_|_|  |___/
##
##

process.filter = cms.Sequence()

## Zskim
process.tagGsf =  cms.EDProducer("CandViewShallowCloneCombiner",
                                 decay = cms.string("PassingWP90 PassingWP90"),
                                 checkCharge = cms.bool(False),
                                 cut   = cms.string("mass > " + str(MASS_CUT_MIN))
                                 )
process.tagGsfCounter = cms.EDFilter("CandViewCountFilter",
                                     src = cms.InputTag("tagGsf"),
                                     minNumber = cms.uint32(1)
                                     )
process.filter *= (process.tagGsf * process.tagGsfCounter)

process.tagGsfSeq = cms.Sequence()

# chiara: questo ancora da rivedere
if (HLTFilter):
    import copy
    from HLTrigger.HLTfilters.hltHighLevel_cfi import *
    process.ZEEHltFilter = copy.deepcopy(hltHighLevel)
    process.ZEEHltFilter.throw = cms.bool(False)
    process.ZEEHltFilter.HLTPaths = ["HLT_Ele*"]
    process.tagGsfSeq *= process.ZEEHltFilter


# chiara: questo ancora da rivedere- secondo me si puo' spegnere tanto rho c'e'
# if (not runFromALCA):
#    process.load('RecoJets.Configuration.RecoPFJets_cff')
#    process.kt6PFJetsForRhoCorrection = process.kt6PFJets.clone(doRhoFastjet = True)
#    process.kt6PFJetsForRhoCorrection.Rho_EtaMax = cms.double(2.5)
#    process.tagGsfSeq *= (process.kt6PFJetsForRhoCorrection) 

process.tagGsfSeq *= (process.ele_sequence * process.filter * process.electronRecalib )

if ( (ECALRecalib) or (not runFromALCA) ):
    process.tagGsfSeq *= ( process.seqALCARECOEcalCalElectronRECO )    

process.zFilterPath = cms.Path( process.tagGsfSeq )

# chiara: quando sistemaranno la def di alcarecooutput rivedere. Per ora (26/11):
# ci aggiungo "keep *_*gedGsfElectron*_*_*" e "keep *_fixedGridRho*_*_*"
# modifico Calibration/EcalAlCaRecoProducers/python/ALCARECOEcalCalIsolElectron_Output_cff.py x salvare i rechits in ES
# rimuovo "keep *_kt6*_rho_*"

process.OutALCARECOEcalCalElectron.outputCommands.extend( [ "keep *_pfMet_*_*", "keep *_offlinePrimaryVerticesWithBS_*_*","keep *_generator_*_*","keep *_*gedGsfElectron*_*_*", "keep *_fixedGridRho*_*_*"] )
if (ECALRecalib):
    process.OutALCARECOEcalCalElectron.outputCommands.extend( [
        "drop recoGsfElectrons_*_*_*",
        "drop recoGsfElectronCores_*_*_*",
        "drop recoSuperClusters_*_*_*",
        "drop recoCaloClusters_*_*_*",
        "drop recoPreshowerClusters_*_*_*",
        "drop recoPreshowerClusterShapes_*_*_*",
        "keep *_electronRecalibSCAssociatorSH_*_*",
        'keep recoSuperClusters_correctedHybridSuperClusters_*_'+processName,
        'keep recoCaloClusters_hybridSuperClusters_*_'+processName,
        'keep recoSuperClusters_hybridSuperClusters_uncleanOnlyHybridSuperClusters_'+processName,
        # Endcap clusters
        'keep recoCaloClusters_multi5x5SuperClusters_multi5x5EndcapBasicClusters_'+processName,
        'keep recoSuperClusters_correctedMulti5x5SuperClustersWithPreshower_*_'+processName,
        # Preshower clusters
        'keep recoPreshowerClusters_multi5x5SuperClustersWithPreshower_*_'+processName,
        'keep recoPreshowerClusterShapes_multi5x5PreshowerClusterShape_*_'+processName ] )
    
process.ALCARECOoutput = cms.OutputModule("PoolOutputModule",
                                          splitLevel = cms.untracked.int32(0),
                                          outputCommands = process.OutALCARECOEcalCalElectron.outputCommands,
                                          fileName = cms.untracked.string('alcaRecoSkim.root'),
                                          SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('zFilterPath')),
                                          dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('ALCARECO')
        )
                                          )                                          
print "OUTPUTCOMMANDS"
print process.ALCARECOoutput.outputCommands
 
process.ALCARECOoutput_step = cms.EndPath(process.ALCARECOoutput)

process.schedule = cms.Schedule(process.zFilterPath,process.ALCARECOoutput_step)
