#include "RecoLocalCalo/EcalRecProducers/plugins/ESDummyRecHitProducer.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <iostream>
#include <cmath>
#include <vector>

#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"

ESDummyRecHitProducer::ESDummyRecHitProducer(const edm::ParameterSet& ps) {

  // std::cout << "ESDummyRecHitProducer::ESDummyRecHitProducer" << std::endl;

  ESRecHitCollection_ = ps.getParameter<edm::InputTag>("ESRecHitCollection");
  ESRecHitToken_ = consumes<ESRecHitCollection>(ESRecHitCollection_);
  ESDummyRecHitCollection_ = ps.getParameter<std::string>("ESDummyRecHitCollection");
  
  produces< ESRecHitCollection >(ESDummyRecHitCollection_);
}

ESDummyRecHitProducer::~ESDummyRecHitProducer() {
  
  // std::cout << "ESDummyRecHitProducer::~ESDummyRecHitProducer" << std::endl;
}

void ESDummyRecHitProducer::produce(edm::Event& evt, const edm::EventSetup& es)
{
  // std::cout << "ESDummyRecHitProducer::produce" << std::endl;
  
  using namespace edm;
  Handle< ESRecHitCollection > pESRecHits;
  
  const ESRecHitCollection*  ESRecHits = 0;
  
  if ( ESRecHitCollection_.label() != "" ) {
    evt.getByToken( ESRecHitToken_, pESRecHits);
    ESRecHits = pESRecHits.product(); // get a ptr to the product
  }
  
  // collection of rechits to put in the event
  std::auto_ptr< ESRecHitCollection > ESDummyRecHits( new ESRecHitCollection );
  
  // loop over rechits to duplicate them
  if (ESRecHits) {
    
    for(ESRecHitCollection::const_iterator it = ESRecHits->begin(); it != ESRecHits->end(); ++it) {
      EcalRecHit aHit( (*it).id(), (*it).energy(), (*it).time() );
      ESDummyRecHits->push_back( aHit );
    }
  }
  
  // put the collection of reconstructed hits in the event   
  LogInfo("EcalDummyRecHitInfo") << "total # ES dummy rechits: " << ESDummyRecHits->size();
  
  evt.put( ESDummyRecHits, ESDummyRecHitCollection_ );
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE( ESDummyRecHitProducer );
