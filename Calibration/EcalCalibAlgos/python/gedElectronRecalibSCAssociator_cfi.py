import FWCore.ParameterSet.Config as cms

gedElectronRecalibSCAssociator = cms.EDProducer("GEDElectronRecalibSuperClusterAssociator",

    gsfElectrons = cms.InputTag('gedGsfElectrons'),
    ebSuperclusters = cms.InputTag('particleFlowSuperClusterECAL:particleFlowSuperClusterECALBarrel'),
    eeSuperclusters = cms.InputTag('particleFlowSuperClusterECAL:particleFlowSuperClusterECALEndcap')
)


