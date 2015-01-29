#ifndef GEDElectronRecalibSuperClusterAssociator_h
#define GEDElectronRecalibSuperClusterAssociator_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronCoreFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronCore.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

#include <string>

class PixelMatchElectronAlgo;

class GEDElectronRecalibSuperClusterAssociator : public edm::EDProducer
{
 public:
  
  explicit GEDElectronRecalibSuperClusterAssociator(const edm::ParameterSet& conf);
  
  virtual ~GEDElectronRecalibSuperClusterAssociator();
  
  virtual void produce(edm::Event& e, const edm::EventSetup& c);
  
 private:

  edm::EDGetTokenT<reco::GsfElectronCollection> electronToken_;     
  edm::EDGetTokenT<reco::SuperClusterCollection> ebScToken_;
  edm::EDGetTokenT<reco::SuperClusterCollection> eeScToken_;
};
#endif
