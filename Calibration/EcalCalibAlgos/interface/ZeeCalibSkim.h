// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"

#include "Math/VectorUtil.h"

class ZeeCalibSkim : public edm::EDFilter {
 public:
  explicit ZeeCalibSkim(const edm::ParameterSet&);
  ~ZeeCalibSkim();
 private:
  virtual void beginJob() {};
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() {};

  // ----------member data ---------------------------
  edm::EDGetTokenT<edm::View<reco::GsfElectron> > electronCollectionToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronIdMapToken_;
  
  float mass_cut_low;
  float mass_cut_high;
 };
