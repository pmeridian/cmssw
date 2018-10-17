/** \file
 *
 */

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/stream/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/EventSetup.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>

#include "RecoMTD/DetLayers/interface/MTDDetLayerGeometry.h"
#include "RecoMTD/Records/interface/MTDRecoGeometryRecord.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"

#include "DataFormats/TrackerRecHit2D/interface/MTDTrackingRecHit.h"

#include "RecoMTD/DetLayers/interface/MTDTrayBarrelLayer.h"
#include "RecoMTD/DetLayers/interface/MTDDetTray.h"
#include "RecoMTD/DetLayers/interface/MTDRingForwardDoubleLayer.h"
#include "RecoMTD/DetLayers/interface/MTDDetRing.h"

#include <DataFormats/MuonDetId/interface/CSCDetId.h>

#include <DataFormats/ForwardDetId/interface/BTLDetId.h>
#include <DataFormats/ForwardDetId/interface/ETLDetId.h>
#include <DataFormats/ForwardDetId/interface/MTDChannelIdentifier.h>
#include "Geometry/CommonTopologies/interface/PixelTopology.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include <sstream>

#include "CLHEP/Random/RandFlat.h"

#include "Geometry/CommonTopologies/interface/Topology.h"

using namespace std;
using namespace edm;

template<class TrackCollection>
class TrackExtenderWithMTDT : public edm::stream::EDProducer<> {  
 public:
  typedef typename TrackCollection:: value_type TrackType;
  typedef edm::View<TrackType> InputCollection;


  TrackExtenderWithMTDT(const ParameterSet& pset); 

  void produce(edm::Event& ev, const edm::EventSetup& es) override final;

  TrackType tryBTLLayers(const TrackType&, 
			 const MTDTrackingDetSetVector&,
			 const MTDDetLayerGeometry*, 
			 const MagneticField* field) const;
  TrackType tryETLLayers(const TrackType&,
			 const MTDTrackingDetSetVector&,
			 const MTDDetLayerGeometry*, 
			 const MagneticField* field) const;

  string dumpLayer(const DetLayer* layer) const;

 private:
  edm::EDGetTokenT<InputCollection> tracksToken_;
  edm::EDGetTokenT<MTDTrackingDetSetVector> hitsToken_;
  std::unique_ptr<MeasurementEstimator> theEstimator;
  edm::ESHandle<TransientTrackBuilder> builder;
};


template<class TrackCollection>  
TrackExtenderWithMTDT<TrackCollection>::TrackExtenderWithMTDT(const ParameterSet& iConfig) {
  float theMaxChi2=25.;
  float theNSigma=3.;
  theEstimator = std::make_unique<Chi2MeasurementEstimator>(theMaxChi2,theNSigma);
  
  tracksToken_ = consumes<InputCollection>(iConfig.getParameter<edm::InputTag>("tracksSrc"));
  hitsToken_ = consumes<MTDTrackingDetSetVector>(iConfig.getParameter<edm::InputTag>("hitsSrc"));

  produces<TrackCollection>();
}

template<class TrackCollection>
void TrackExtenderWithMTDT<TrackCollection>::produce( edm::Event& ev,
						      const edm::EventSetup& es ) {
  
  edm::ESHandle<MTDDetLayerGeometry> geo;
  es.get<MTDRecoGeometryRecord>().get(geo);

  edm::ESHandle<MagneticField> magfield;
  es.get<IdealMagneticFieldRecord>().get(magfield);  
    
  es.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);

  // Some printouts

  auto output = std::make_unique<TrackCollection>();

  cout << "*** allBTLLayers(): " << geo->allBTLLayers().size() << endl;
  for (auto dl = geo->allBTLLayers().begin();
       dl != geo->allBTLLayers().end(); ++dl) {
    cout << "  " << (int) (dl-geo->allBTLLayers().begin()) << " " << dumpLayer(*dl);
  }
  cout << endl << endl;
  
  cout << "*** allETLLayers(): " << geo->allETLLayers().size() << endl;
  for (auto dl = geo->allETLLayers().begin();
       dl != geo->allETLLayers().end(); ++dl) {
    cout << "  " << (int) (dl-geo->allETLLayers().begin()) << " " << dumpLayer(*dl);
  }
  cout << endl << endl;
  
  cout << "*** allLayers(): " << geo->allLayers().size() << endl;
  for (auto dl = geo->allLayers().begin();
       dl != geo->allLayers().end(); ++dl) {
    cout << "  " << (int) (dl-geo->allLayers().begin()) << " " << dumpLayer(*dl);
  }
  cout << endl << endl;

  edm::Handle<InputCollection> tracksH;  
  ev.getByToken(tracksToken_,tracksH);
  const auto& tracks = *tracksH;

  edm::Handle<MTDTrackingDetSetVector> hitsH;  
  ev.getByToken(hitsToken_,hitsH);
  const auto& hits = *hitsH;

  for( const auto& track : tracks ) {  
    std::cout << "NEW TRACK" << std::endl;
    tryBTLLayers(track,hits,geo.product(),magfield.product());
    tryETLLayers(track,hits,geo.product(),magfield.product());
  }

  ev.put(std::move(output));
}

template<class TrackCollection>
typename TrackExtenderWithMTDT<TrackCollection>::TrackType
TrackExtenderWithMTDT<TrackCollection>::tryBTLLayers(const TrackType& track,
						     const MTDTrackingDetSetVector& hits,
						     const MTDDetLayerGeometry* geo,
						     const MagneticField* field) const {
  auto output = track;
  const vector<const DetLayer*>& layers = geo->allBTLLayers();

  auto cmp = [](const unsigned one, const unsigned two) -> bool { return one < two; };

  auto tTrack = builder->build(track);

  for (auto ilay = layers.begin(); ilay!=layers.end(); ++ilay) {
    const MTDTrayBarrelLayer* layer = (const MTDTrayBarrelLayer*) (*ilay);
    
    // get the outermost trajectory point on the track    
    TrajectoryStateOnSurface tsos = tTrack.outermostMeasurementState();
    cout << "testBTLLayers: at " << tsos.globalPosition()
	 << " R=" << tsos.globalPosition().perp()
	 << " phi=" << tsos.globalPosition().phi()
	 << " Z=" << tsos.globalPosition().z()
	 << " p = " << tsos.globalMomentum()
	 << endl;
    
    SteppingHelixPropagator prop(field,anyDirection);

    pair<bool, TrajectoryStateOnSurface> comp = layer->compatible(tsos,prop,*theEstimator);
    if( comp.first ) {
      cout << "is compatible: " << comp.first
	   << " at: R=" << comp.second.globalPosition().perp()
	   << " phi=" << comp.second.globalPosition().phi()
	   << " Z=" <<  comp.second.globalPosition().z()
	   << endl;

    
      vector<DetLayer::DetWithState> compDets = layer->compatibleDets(tsos,prop,*theEstimator);
      if (compDets.size()) {
	for( const auto& detWithState : compDets ) {
	  
	  auto prop_pos = detWithState.second.globalPosition();
	  auto local_pos = detWithState.first->toLocal(prop_pos);
	  auto meas_point = detWithState.first->topology().measurementPosition(local_pos);
	  auto pixel = static_cast<const PixelTopology&>(detWithState.first->topology()).pixel(local_pos);
	  auto theChannel =  detWithState.first->topology().channel(local_pos);
	  auto pixFromCh = MTDChannelIdentifier::channelToPixel(theChannel);
	  
	  
	  cout << "compatibleDets: " << compDets.size() << endl
	       << "  final state pos: " << prop_pos << endl 
	       << "  local state pos: " << local_pos << endl
	       << "  measurement pos: " << meas_point << endl
	       << " pixel : (" << pixel.first << ',' << pixel.second << ')' << endl
	       << "       : (" << pixFromCh.first << ',' << pixFromCh.second << ')' << endl
	       << "  channel : " << theChannel << endl
	       << "  det         pos: " << detWithState.first->position()
	       << " id: " << std::hex << BTLDetId(detWithState.first->geographicalId().rawId()).rawId() << std::dec<< endl 
	       << "  distance " << (tsos.globalPosition()-detWithState.first->position()).mag()
	       << endl
	       << endl; 

	  auto range = hits.equal_range(detWithState.first->geographicalId(),cmp);
	  for( auto detitr = range.first; detitr != range.second; ++detitr ) {
	    std::cout << "\tdet with hits: " << std::distance(range.first,detitr) << ' ' << detitr->size() << std::endl;
	    for( auto itr = detitr->begin(); itr != detitr->end(); ++itr ) {
	      std::cout << itr->localPosition() << ' ' << itr->localPositionError() << std::endl;
	    }
	  }	  
	}     
	
      } else {
	cout << " ERROR : no compatible BTL det found" << endl;
      }    
    }
  }
  return output;
}

template<class TrackCollection>
typename TrackExtenderWithMTDT<TrackCollection>::TrackType
TrackExtenderWithMTDT<TrackCollection>::tryETLLayers(const TrackType& track, 
						     const MTDTrackingDetSetVector& hits,
						     const MTDDetLayerGeometry* geo,
						     const MagneticField* field) const {
  auto output = track;
  const vector<const DetLayer*>& layers = geo->allETLLayers();

  auto cmp = [](const unsigned one, const unsigned two) -> bool { return one < two; };

  auto tTrack = builder->build(track);

  for (auto ilay = layers.begin(); ilay!=layers.end(); ++ilay) {
    const MTDRingForwardDoubleLayer* layer = (const MTDRingForwardDoubleLayer*) (*ilay);
    const BoundDisk& disk = layer->specificSurface();
    const double diskZ = disk.position().z();

    // get the outermost trajectory point on the track    
    TrajectoryStateOnSurface tsos = tTrack.outermostMeasurementState();
    cout << "testETLLayers: at " << tsos.globalPosition()
	 << " R=" << tsos.globalPosition().perp()
	 << " phi=" << tsos.globalPosition().phi()
	 << " Z=" << tsos.globalPosition().z()
	 << " p = " << tsos.globalMomentum()
	 << endl;

    if( tsos.globalPosition().z() * diskZ < 0 ) continue; // only propagate to the disk that's on the same side

    SteppingHelixPropagator prop(field,anyDirection);

    pair<bool, TrajectoryStateOnSurface> comp = layer->compatible(tsos,prop,*theEstimator);
    if( comp.first ) {
      cout << "is compatible: " << comp.first
	   << " at: R=" << comp.second.globalPosition().perp()
	   << " phi=" << comp.second.globalPosition().phi()
	   << " Z=" <<  comp.second.globalPosition().z()
	   << endl;
      auto gp = comp.second.globalPosition();
            
      // if we're compatible with the disc, try to find modules

      vector<DetLayer::DetWithState> compDets = layer->compatibleDets(tsos,prop,*theEstimator);
      const bool hasDets = compDets.size();
      if( hasDets ) {
	for( const auto& detWithState : compDets ) {
	  auto prop_pos = detWithState.second.globalPosition();
	  auto local_pos = detWithState.first->toLocal(prop_pos);
	  auto meas_point = detWithState.first->topology().measurementPosition(local_pos);
	  auto pixel = static_cast<const PixelTopology&>(detWithState.first->topology()).pixel(local_pos);
	  auto theChannel =  detWithState.first->topology().channel(local_pos);
	  auto pixFromCh = MTDChannelIdentifier::channelToPixel(theChannel);
	  
	  cout << "compatibleDets: " << compDets.size() << endl
	       << "  final state pos: " << prop_pos << endl 
	       << "  local state pos: " <<  local_pos << endl
	       << "  measurement pos: " << meas_point << endl
	       << "  channel: " << theChannel << endl
	       << " pixel : (" << pixel.first << ',' << pixel.second << ')' << endl
	       << "       : (" << pixFromCh.first << ',' << pixFromCh.second << ')' << endl
	       << "  det         pos: " << detWithState.first->position()
	       << " id: " << std::hex << ETLDetId(detWithState.first->geographicalId().rawId()).rawId() << std::dec << endl 
	       << "  distance " << (tsos.globalPosition()-detWithState.first->position()).mag()
	       << endl
	       << endl;
	  
	  auto range = hits.equal_range(detWithState.first->geographicalId(),cmp);
	  for( auto detitr = range.first; detitr != range.second; ++detitr ) {
	    std::cout << "\tdet with hits: " << std::distance(range.first,detitr) << ' ' << detitr->size() << std::endl;
	    for( auto itr = detitr->begin(); itr != detitr->end(); ++itr ) {
	      std::cout << itr->localPosition() << ' ' << itr->localPositionError() << std::endl;
	    }
	  }
	  
	}
      } else {
	if(layer->isCrack(gp))
	  {
	    cout << " CSC crack found " << std::endl;
	  }
	else
	  {
	    cout << " ERROR : no compatible det found in ETL"
		 << " at: R=" << gp.perp()
		 << " phi= " << gp.phi().degrees()
		 << " Z= " << gp.z() << std::endl;
	  }

      }
    }
  }
  return output;
}

template<class TrackCollection>
string TrackExtenderWithMTDT<TrackCollection>::dumpLayer(const DetLayer* layer) const {
  stringstream output;
  
  const BoundSurface* sur=0;
  const BoundCylinder* bc=0;
  const BoundDisk* bd=0;

  sur = &(layer->surface());
  if ( (bc = dynamic_cast<const BoundCylinder*>(sur)) ) {
    output << "  Cylinder of radius: " << bc->radius() << endl;
  }
  else if ( (bd = dynamic_cast<const BoundDisk*>(sur)) ) {
    output << "  Disk at: " <<  bd->position().z() << endl;
  }
  return output.str();
}

//define this as a plug-in
#include <FWCore/Framework/interface/MakerMacros.h>
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
typedef TrackExtenderWithMTDT<reco::TrackCollection> TrackExtenderWithMTD;
typedef TrackExtenderWithMTDT<reco::GsfTrackCollection> GSFTrackExtenderWithMTD;

DEFINE_FWK_MODULE(TrackExtenderWithMTD);
DEFINE_FWK_MODULE(GSFTrackExtenderWithMTD);
