import FWCore.ParameterSet.Config as cms

# From which kind of dataset are you starting from?
MC = True

#electron cuts
ELECTRON_ET_CUT_MIN = 20.0
ELECTRON_CUTS = "(abs(superCluster.eta)<2.5) && (ecalEnergy*sin(superClusterPosition.theta)>" + str(ELECTRON_ET_CUT_MIN) + ")"

MASS_CUT_MIN = 60.

process = cms.Process("ZEESKIM")

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(200)
)

readFiles = cms.untracked.vstring()
readFiles.extend( [
        "/store/user/crovelli/testCalibZee/fileAOD.root"   # MC DY file from phys14
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


process.selectedElectrons = cms.EDFilter("GsfElectronRefSelector",
                                 src = cms.InputTag( 'gedGsfElectrons' ),
                                 cut = cms.string( ELECTRON_CUTS )
                             )

process.ele_sequence = cms.Sequence(
    process.selectedElectrons
    )


process.filter = cms.Sequence()

process.tagGsf =  cms.EDProducer("CandViewShallowCloneCombiner",
                                 decay = cms.string("selectedElectrons selectedElectrons"),
                                 checkCharge = cms.bool(False),
                                 cut   = cms.string("mass > " + str(MASS_CUT_MIN))
                                 )

process.tagGsfCounter = cms.EDFilter("CandViewCountFilter",
                                     src = cms.InputTag("tagGsf"),
                                     minNumber = cms.uint32(1)
                                     )

process.filter *= (process.tagGsf * process.tagGsfCounter)

process.tagGsfSeq = cms.Sequence()

process.tagGsfSeq *= (process.ele_sequence * process.filter)

process.zFilterPath = cms.Path( process.tagGsfSeq )

process.AODEventContent.outputCommands.extend( [
        "drop *_*_*_*",
        "keep *_*gedGsfElectron*_*_*",
        "keep *_*reducedEcalRecHits*_*_*",
        "keep *_*TriggerResults*_*_HLT",
        "keep *_offlinePrimaryVertices_*_*",    # we do not save the collection with BS for the moment
        "keep *_*conversions*_*_*",
        "keep *_*offlineBeamSpot*_*_*",
        "keep *_*particleFlowEGamma*_*_*",
        "keep *_*particleFlowSuperClusterECAL*_*_*",
        "keep *_*electronGsfTracks*_*_*",
        "keep *recoGenParticles_*genParticles*_*_*"] )

process.ZeeSkimOutput = cms.OutputModule("PoolOutputModule",
                                         splitLevel = cms.untracked.int32(0),
                                         outputCommands = process.AODEventContent.outputCommands,
                                         fileName = cms.untracked.string('zeeSkimChiara.root'),
                                         SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('zFilterPath')),
                                         dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('ALCARECO')
        )
                                          )                                          
#print "OUTPUTCOMMANDS"
#print process.AODEventContent.outputCommands
 
process.ZeeSkimOutput_step = cms.EndPath(process.ZeeSkimOutput)

process.schedule = cms.Schedule(process.zFilterPath,process.ZeeSkimOutput_step)
