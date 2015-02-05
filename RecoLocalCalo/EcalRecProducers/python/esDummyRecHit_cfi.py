import FWCore.ParameterSet.Config as cms

# dummy ES rechit producer
ecalPreshowerRecHit = cms.EDProducer("ESDummyRecHitProducer",
    ESRecHitCollection = cms.InputTag("ecalPreshowerRecHit","EcalRecHitsES"),
    ESDummyRecHitCollection = cms.string('EcalRecHitsES')
)
