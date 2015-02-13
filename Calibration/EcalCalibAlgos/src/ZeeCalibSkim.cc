#include "Calibration/EcalCalibAlgos/interface/ZeeCalibSkim.h"

#include <iostream>
#include <utility>

ZeeCalibSkim::ZeeCalibSkim(const edm::ParameterSet& iConfig):
  electronCollectionToken_(consumes<edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("electrons"))),
  electronIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronIdMap"))),
  mass_cut_low(static_cast<float>(iConfig.getUntrackedParameter<double>("mass_cut_low",60.))),
  mass_cut_high(static_cast<float>(iConfig.getUntrackedParameter<double>("mass_cut_high",9999.))),
  requireOppositeCharge(static_cast<bool>(iConfig.getUntrackedParameter<bool>("requireOppositeCharge",false)))
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
  
  std::vector<const reco::GsfElectron*> selElectrons;
  for( const auto& el : ele_refs ) 
    {
      bool isPass;
      isPass = (*id_decisions)[el];
      if (isPass)
	selElectrons.push_back(&(*el));
    }
  
  if (selElectrons.size()<2)
    return false;
  
  std::vector<const reco::GsfElectron*>::const_iterator ibegin=selElectrons.begin(),
    iend = selElectrons.end(), iele = ibegin, jele = iele + 1;
  for ( ; iele != iend - 1; ++iele ) {
    for ( ; jele != iend; ++jele ) {

      if (requireOppositeCharge && (*iele)->charge()+(*jele)->charge() != -1)
	continue;

      //di-electron object
      Candidate::LorentzVector diEle=(*iele)->p4()+(*jele)->p4();

      if (diEle.M()>=mass_cut_low &&
	  diEle.M()<=mass_cut_high
	  )
	  return true;
    }
  }
  return false;
}
