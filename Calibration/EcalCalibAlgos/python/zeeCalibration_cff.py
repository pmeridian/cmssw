import FWCore.ParameterSet.Config as cms

looper = cms.Looper("ZeeCalibration",
    maxLoops = cms.untracked.uint32(10),
    wantEtaCorrection = cms.untracked.bool(False),
    minInvMassCut = cms.untracked.double(70),
    maxInvMassCut = cms.untracked.double(110),
    electronSelection = cms.untracked.int32(-1),
    histoFile = cms.string('myHistograms.root'),
    zeeFile = cms.string('myHistograms.root'),
    initialMiscalibrationBarrel = cms.untracked.string(''),
    initialMiscalibrationEndcap = cms.untracked.string(''),
    etaBins = cms.untracked.uint32(10),
    etBins = cms.untracked.uint32(10),
    etaMin = cms.untracked.double(0.),
    etMin = cms.untracked.double(0.),
    etaMax = cms.untracked.double(3.),
    etMax = cms.untracked.double(100.),
    ZCalib_CalibType = cms.untracked.string('RING'),
    ZCalib_InvMass = cms.untracked.string('SCMass'),                
    #
    rechitProducer = cms.string('recalibRechit'),
    rechitCollection = cms.string('EcalRecHitsEB'),
    erechitProducer = cms.string('recalibRechit'),
    erechitCollection = cms.string('EcalRecHitsEE'),
    ebSuperclusters = cms.InputTag("particleFlowSuperClusterECAL","recalibParticleFlowSuperClusterECALBarrel","ZEECALIB"),
    eeSuperclusters = cms.InputTag("particleFlowSuperClusterECAL","recalibParticleFlowSuperClusterECALEndcapWithPreshower","ZEECALIB"),
    electrons = cms.InputTag("gedGsfElectrons","","RECO"),
    vertices = cms.InputTag("offlinePrimaryVertices"),
    conversions = cms.InputTag("conversions"),
    beamspot = cms.InputTag("offlineBeamSpot"),
    mcProducer = cms.untracked.string('genParticles'),
    HLTriggerResults = cms.InputTag("TriggerResults","","HLT")
)



