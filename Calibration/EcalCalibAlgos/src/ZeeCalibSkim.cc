#include "Calibration/EcalCalibAlgos/interface/ZeeCalibSkim.h"

#include <iostream>

ZeeCalibSkim::ZeeCalibSkim(const edm::ParameterSet& iConfig):
  electronCollectionToken_(consumes<edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("electrons"))),
  electronIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronIdMap")))
{
}

ZeeCalibSkim::~ZeeCalibSkim()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}


// ------------ method called for each event ------------
bool
ZeeCalibSkim::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  using namespace reco;

  // Get the electron collection
  edm::Handle<edm::View<reco::GsfElectron> > collection;
  iEvent.getByToken(electronCollectionToken_, collection);

  // Get the electron ID data from the event stream.
  // Note: this implies that the VID ID modules have been run upstream.
  edm::Handle<edm::ValueMap<bool> > id_decisions;
  iEvent.getByToken(electronIdMapToken_,id_decisions);
  
  // Loop over electrons
  const auto& ele_refs = collection->refVector();
  for( const auto& el : ele_refs ) {
    bool isPass;
    isPass = (*id_decisions)[el];
    std::cout << isPass << std::endl;
  }
  return true;
}
