import FWCore.ParameterSet.Config as cms

from RecoLocalFastTime.FTLRecProducers.ftlUncalibratedRecHits_cfi import ftlUncalibratedRecHits
from RecoLocalFastTime.FTLRecProducers.ftlRecHits_cfi import ftlRecHits

fastTimingLocalReco = cms.Sequence(ftlUncalibratedRecHits*ftlRecHits)

from RecoLocalFastTime.FTLRecProducers.mtdUncalibratedRecHits_cfi import mtdUncalibratedRecHits
from RecoLocalFastTime.FTLRecProducers.mtdRecHits_cfi import mtdRecHits

_phase2_timing_layer_tile_fastTimingLocalReco = cms.Sequence(mtdUncalibratedRecHits*mtdRecHits*mtdClusters*mtdTrackingRecHits)
from Configuration.Eras.Modifier_phase2_timing_layer_tile_cff import phase2_timing_layer_tile
phase2_timing_layer_tile.toReplaceWith(fastTimingLocalReco, _phase2_timing_layer_tile_fastTimingLocalReco)
phase2_timing_layer_tile.toModify(mtdRecHits, barrelUncalibratedRecHits = cms.InputTag('mtdUncalibratedRecHits:FTLBarrel'),
                                  endcapUncalibratedRecHits = cms.InputTag('mtdUncalibratedRecHits:FTLEndcap') )

_phase2_timing_layer_bar_fastTimingLocalReco = cms.Sequence(mtdUncalibratedRecHits*mtdRecHits*mtdClusters*mtdTrackingRecHits)
from Configuration.Eras.Modifier_phase2_timing_layer_bar_cff import phase2_timing_layer_bar
phase2_timing_layer_bar.toReplaceWith(fastTimingLocalReco, _phase2_timing_layer_bar_fastTimingLocalReco)
phase2_timing_layer_bar.toModify(mtdRecHits, barrelUncalibratedRecHits = cms.InputTag('mtdUncalibratedRecHits:FTLBarrel'),
                                 endcapUncalibratedRecHits = cms.InputTag('mtdUncalibratedRecHits:FTLEndcap') )


