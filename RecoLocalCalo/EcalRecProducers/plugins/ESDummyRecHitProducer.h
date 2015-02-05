#ifndef RecoLocalCalo_EcalRecProducers_ESDummyRecHitProducer_HH
#define RecoLocalCalo_EcalRecProducers_ESDummyRecHitProducer_HH

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

class ESDummyRecHitProducer : public edm::EDProducer {

        public:
                explicit ESDummyRecHitProducer(const edm::ParameterSet& ps);
                ~ESDummyRecHitProducer();
                virtual void produce(edm::Event& evt, const edm::EventSetup& es);

        private:

		edm::InputTag ESRecHitCollection_;
		edm::EDGetTokenT<ESRecHitCollection> ESRecHitToken_;
                std::string ESDummyRecHitCollection_; 
};
#endif
