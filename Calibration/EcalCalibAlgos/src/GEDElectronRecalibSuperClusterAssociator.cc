// -*- C++ -*-
//
//

// user include files
#include "Calibration/EcalCalibAlgos/interface/GEDElectronRecalibSuperClusterAssociator.h"

#include <iostream>

#define DEBUG

using namespace reco;
using namespace edm;
using namespace std;

GEDElectronRecalibSuperClusterAssociator::GEDElectronRecalibSuperClusterAssociator(const edm::ParameterSet& iConfig)  
{
#ifdef DEBUG
  std::cout<< "GEDElectronRecalibSuperClusterAssociator::GEDElectronRecalibSuperClusterAssociator" << std::endl;
#endif

  electronToken_ = consumes<GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("gsfElectrons"));
  ebScToken_ = consumes<SuperClusterCollection>(iConfig.getParameter<edm::InputTag>("ebSuperclusters"));
  eeScToken_ = consumes<SuperClusterCollection>(iConfig.getParameter<edm::InputTag>("eeSuperclusters"));

  // register your products
  produces<GsfElectronCollection>();
  produces<GsfElectronCoreCollection>() ;
  
#ifdef DEBUG
  std::cout<< "GEDElectronRecalibSuperClusterAssociator::GEDElectronRecalibSuperClusterAssociator::end" << std::endl;
#endif
}

GEDElectronRecalibSuperClusterAssociator::~GEDElectronRecalibSuperClusterAssociator() { }

void GEDElectronRecalibSuperClusterAssociator::produce(edm::Event& e, const edm::EventSetup& iSetup)
{
#ifdef DEBUG
  std::cout<< "GEDElectronRecalibSuperClusterAssociator::produce" << std::endl;
#endif

  // Create the output collections
  std::auto_ptr<GsfElectronCollection> pOutEle(new GsfElectronCollection);
  std::auto_ptr<GsfElectronCoreCollection> pOutEleCore(new GsfElectronCoreCollection);


  // Get SuperClusters in EB
  Handle<reco::SuperClusterCollection> ebScHandle;
  e.getByToken(ebScToken_, ebScHandle); 
  const reco::SuperClusterCollection* ebScCollection = ebScHandle.product();

#ifdef DEBUG
  std::cout<<"GEDElectronRecalibSuperClusterAssociator::EB scCollection->size()" << ebScCollection->size()<<std::endl;
#endif
 
  // Get SuperClusters in EE
  Handle<reco::SuperClusterCollection> eeScHandle;
  e.getByToken(eeScToken_, eeScHandle); 
  const reco::SuperClusterCollection* eeScCollection = eeScHandle.product();

#ifdef DEBUG
  std::cout<<"GEDElectronRecalibSuperClusterAssociator::EE scCollection->size()" << eeScCollection->size() << std::endl;
#endif

  // Get Electrons
  edm::Handle<reco::GsfElectronCollection> eleHandle;
  e.getByToken(electronToken_, eleHandle);
  const reco::GsfElectronCollection* electronCollection = eleHandle.product();

#ifdef DEBUG
  std::cout<<"GEDElectronRecalibSuperClusterAssociator::Electron Collection->size()" << electronCollection->size() << std::endl;
#endif

  GsfElectronCoreRefProd rEleCore=e.getRefBeforePut<GsfElectronCoreCollection>();
  edm::Ref<GsfElectronCoreCollection>::key_type idxEleCore = 0;

  for(reco::GsfElectronCollection::const_iterator eleIt = electronCollection->begin(); eleIt != electronCollection->end(); eleIt++) {

    float DeltaRMineleSCbarrel(0.15);
    float DeltaRMineleSCendcap(0.15);
    const reco::SuperCluster* nearestSCbarrel=0;
    const reco::SuperCluster* nearestSCendcap=0;
    int iSC=0;

    // loop on EB superClusters
    int iscRefEB=-1;
    for(reco::SuperClusterCollection::const_iterator scEbIt = ebScCollection->begin(); scEbIt != ebScCollection->end(); scEbIt++) {
	
#ifdef DEBUG
      std::cout << "GEDElectronRecalibSuperClusterAssociator::sc in EB" << endl; 
      std::cout << scEbIt->energy() << " " << scEbIt->eta() << " " << scEbIt->phi() << " " << eleIt->eta() << " " << eleIt->phi() << std::endl;
#endif
      double DeltaReleSC = sqrt ( pow( eleIt->eta() - scEbIt->eta(),2) + pow(eleIt->phi() - scEbIt->phi(),2));
      if(DeltaReleSC<DeltaRMineleSCbarrel) {
	DeltaRMineleSCbarrel = DeltaReleSC;
	nearestSCbarrel = &*scEbIt;
	iscRefEB = iSC;
      }
      iSC++;
    }
    iSC = 0;

    // loop on EE superClusters
    int iscRefEE=-1;
    for(reco::SuperClusterCollection::const_iterator scEeIt = eeScCollection->begin(); scEeIt != eeScCollection->end(); scEeIt++){
      
#ifdef DEBUG
      std::cout << "GEDElectronRecalibSuperClusterAssociator::sc in EE" << endl; 
      std::cout << "EE " << scEeIt->energy() << " " << scEeIt->eta() << " " << scEeIt->phi() << " " << eleIt->eta() << " " << eleIt->phi() << std::endl;
#endif
      double DeltaReleSC = sqrt ( pow( eleIt->eta() - scEeIt->eta(),2) + pow(eleIt->phi() - scEeIt->phi(),2));
      if(DeltaReleSC<DeltaRMineleSCendcap) {
	DeltaRMineleSCendcap = DeltaReleSC;
	nearestSCendcap = &*scEeIt;
	iscRefEE = iSC;
      }
      iSC++;
    }
  
#ifdef DEBUG
    std::cout << "+++++++++++" << std::endl;
    std::cout << &(*eleIt->gsfTrack()) << std::endl;
    std::cout << &(*eleIt->superCluster()) << std::endl;
    std::cout << "+++++++++++" << std::endl;
#endif
  
    if(eleIt->isEB() && nearestSCbarrel){
      
#ifdef DEBUG
      std::cout << "Starting Association is with EB superCluster "<< std::endl;
#endif 
      
      reco::GsfElectronCore newEleCore(*(eleIt->core()));
      newEleCore.setGsfTrack(eleIt->gsfTrack());                 // gsf track
      std::cout << "newEleCore.setGsfTrack done" << std::endl;
      newEleCore.setSuperCluster(eleIt->superCluster());         // refined supercluster
      std::cout << "newEleCore.setSuperCluster done" << std::endl;
      reco::SuperClusterRef scRef(reco::SuperClusterRef(ebScHandle, iscRefEB));    
      newEleCore.setParentSuperCluster(scRef);                   // mustache 
      std::cout << "newEleCore.setParentSuperCluster done" << std::endl;
      std::cout << "+++++++  ==> " <<  &(*scRef) << " <== +++++++" << endl;

      reco::GsfElectronCoreRef newEleCoreRef(reco::GsfElectronCoreRef(rEleCore, idxEleCore ++));
      // reco::GsfElectronCoreRef newEleCoreRef(reco::GsfElectronCoreRef(newEleCore));
      std::cout << "creato il newEleCoreRef" << std::endl;
      std::cout << "+++++++  +++ ==> " <<  &(*newEleCoreRef) << " <== +++ +++++++" << endl;

      pOutEleCore->push_back(newEleCore);
      std::cout << "fatto il push back del core" << std::endl;

      // chiara: qui e' dove crasha
      reco::GsfElectron newEle(*eleIt,newEleCoreRef);
      // reco::GsfElectron newEle(*eleIt,newEleCore);
      std::cout << "ok new Ele" << endl;
      // chiara: qui e' dove crasha  

      // pOutEleCore->push_back(newEleCore);
      // std::cout << "fatto il push back del core" << std::endl;

      // chiara: da vedere se va bene cosi' (specie per le correzioni)
      newEle.setP4( eleIt->p4()*(nearestSCbarrel->energy()/eleIt->parentSuperCluster()->energy()) );
      // newEle.setCorrectedEcalEnergy(eleIt->p4().energy()*(nearestSCbarrel->energy()/eleIt->ecalEnergy()),eleIt->ecalEnergyError()*(nearestSCbarrel->energy()/eleIt->ecalEnergy()));             

      pOutEle->push_back(newEle);
      
#ifdef DEBUG
      std::cout << "Association is with EB superCluster "<< std::endl;
#endif 
    }

    if(!(eleIt->isEB()) && nearestSCendcap) {
      
#ifdef DEBUG
      std::cout << "Starting Association is with EE superCluster "<< std::endl;
#endif
      
      reco::GsfElectronCore newEleCore(*(eleIt->core()));
      newEleCore.setGsfTrack(eleIt->gsfTrack());                 // gsf track
      newEleCore.setSuperCluster(eleIt->superCluster());         // refined supercluster
      reco::SuperClusterRef scRef(reco::SuperClusterRef(eeScHandle, iscRefEE));  
      newEleCore.setParentSuperCluster(scRef);                   // mustache 
      reco::GsfElectronCoreRef newEleCoreRef(reco::GsfElectronCoreRef(rEleCore, idxEleCore ++));
      pOutEleCore->push_back(newEleCore);
      reco::GsfElectron newEle(*eleIt,newEleCoreRef);

      // chiara: da vedere se va bene cosi' (specie per le correzioni)
      newEle.setP4( eleIt->p4()*(nearestSCendcap->energy()/eleIt->parentSuperCluster()->energy()) );
      // newEle.setCorrectedEcalEnergy(eleIt->p4().energy()*(nearestSCendcap->energy()/eleIt->ecalEnergy()),eleIt->ecalEnergyError()*(nearestSCendcap->energy()/eleIt->ecalEnergy()));
      // std::cout << "FROM REF " << newEle.superCluster().key() << std::endl;
      pOutEle->push_back(newEle);
      
#ifdef DEBUG
      std::cout << "Association is with EE superCluster "<< std::endl;
#endif
    }
  }
  
  
#ifdef DEBUG
  std::cout << "Filled new electrons "     << pOutEle->size()     << std::endl;
  std::cout << "Filled new electronsCore " << pOutEleCore->size() << std::endl;
#endif

  // put result into the Event
  e.put(pOutEle);
  e.put(pOutEleCore);
} 

