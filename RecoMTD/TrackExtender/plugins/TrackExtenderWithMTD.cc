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
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
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

#include "TrackingTools/PatternTools/interface/Trajectory.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"

#include "RecoMTD/TransientTrackingRecHit/interface/MTDTransientTrackingRecHitBuilder.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"

#include "TrackingTools/PatternTools/interface/TSCBLBuilderWithPropagator.h"

#include "RecoTracker/TransientTrackingRecHit/interface/Traj2TrackHits.h"
#include "TrackingTools/TrackRefitter/interface/TrackTransformer.h"

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

  TransientTrackingRecHit::ConstRecHitContainer tryBTLLayers(const TrackType&, 
							     const MTDTrackingDetSetVector&,
							     const MTDDetLayerGeometry*, 
							     const MagneticField* field) const;

  TransientTrackingRecHit::ConstRecHitContainer tryETLLayers(const TrackType&,
							     const MTDTrackingDetSetVector&,
							     const MTDDetLayerGeometry*, 
							     const MagneticField* field) const;
  
  RefitDirection::GeometricalDirection
  checkRecHitsOrdering(TransientTrackingRecHit::ConstRecHitContainer const & recHits) const {
    
    if (!recHits.empty()){
      GlobalPoint first = gtg->idToDet(recHits.front()->geographicalId())->position();
      GlobalPoint last = gtg->idToDet(recHits.back()->geographicalId())->position();
      
      // maybe perp2?
      auto rFirst = first.mag2();
      auto rLast  = last.mag2();
      if(rFirst < rLast) return RefitDirection::insideOut;
      if(rFirst > rLast) return RefitDirection::outsideIn;
    }
    LogDebug("Reco|TrackingTools|TrackTransformer") << "Impossible to determine the rechits order" <<endl;
    return RefitDirection::undetermined;
  }

  reco::Track buildTrack(const reco::Track&, const Trajectory&, const reco::BeamSpot&, const MagneticField* field, bool hasMTD) const;
  reco::TrackExtra buildTrackExtra(const Trajectory& trajectory) const;

  string dumpLayer(const DetLayer* layer) const;

 private:
  edm::EDGetTokenT<InputCollection> tracksToken_;
  edm::EDGetTokenT<MTDTrackingDetSetVector> hitsToken_;
  edm::EDGetTokenT<reco::BeamSpot> bsToken_;
  const bool updateTraj_, updateExtra_, updatePattern_;
  std::unique_ptr<MeasurementEstimator> theEstimator;
  std::unique_ptr<TrackTransformer> theTransformer;
  edm::ESHandle<TransientTrackBuilder> builder;
  edm::ESHandle<TransientTrackingRecHitBuilder> hitbuilder;
  edm::ESHandle<GlobalTrackingGeometry> gtg;
};


template<class TrackCollection>  
TrackExtenderWithMTDT<TrackCollection>::TrackExtenderWithMTDT(const ParameterSet& iConfig) :
  updateTraj_(iConfig.getParameter<bool>("updateTrackTrajectory")),
  updateExtra_(iConfig.getParameter<bool>("updateTrackExtra")),
  updatePattern_(iConfig.getParameter<bool>("updateTrackHitPattern")) {
  float theMaxChi2=25.;
  float theNSigma=3.;
  theEstimator = std::make_unique<Chi2MeasurementEstimator>(theMaxChi2,theNSigma);
  
  theTransformer = std::make_unique<TrackTransformer>(iConfig.getParameterSet("TrackTransformer"));

  tracksToken_ = consumes<InputCollection>(iConfig.getParameter<edm::InputTag>("tracksSrc"));
  hitsToken_ = consumes<MTDTrackingDetSetVector>(iConfig.getParameter<edm::InputTag>("hitsSrc"));
  bsToken_ = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotSrc"));

  produces<edm::OwnVector<TrackingRecHit>>();
  produces<reco::TrackExtraCollection>();
  produces<TrackCollection>();
}

template<class TrackCollection>
void TrackExtenderWithMTDT<TrackCollection>::produce( edm::Event& ev,
						      const edm::EventSetup& es ) {
  
  //this produces pieces of the track extra
  Traj2TrackHits t2t;

  theTransformer->setServices(es);

  TrackingRecHitRefProd hitsRefProd = ev.getRefBeforePut<TrackingRecHitCollection>();
  reco::TrackExtraRefProd extrasRefProd = ev.getRefBeforePut<reco::TrackExtraCollection>();
  
  es.get<GlobalTrackingGeometryRecord>().get(gtg);

  edm::ESHandle<MTDDetLayerGeometry> geo;
  es.get<MTDRecoGeometryRecord>().get(geo);

  edm::ESHandle<MagneticField> magfield;
  es.get<IdealMagneticFieldRecord>().get(magfield);  
    
  es.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);

  es.get<TransientRecHitRecord>().get("MTDRecHitBuilder",hitbuilder);

  edm::ESHandle<TrackerTopology> httopo;
  es.get<TrackerTopologyRcd>().get(httopo);
  const TrackerTopology& ttopo = *httopo;

  // Some printouts

  auto output  = std::make_unique<TrackCollection>();
  auto extras  = std::make_unique<reco::TrackExtraCollection>();
  auto outhits = std::make_unique<edm::OwnVector<TrackingRecHit>>();

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

  edm::Handle<reco::BeamSpot> bsH;  
  ev.getByToken(bsToken_,bsH);
  const auto& bs = *bsH;

  for( const auto& track : tracks ) {  
    std::cout << "NEW TRACK" << std::endl;
    reco::TransientTrack ttrack(track,magfield.product(),gtg);
    auto trajs = theTransformer->transform(track);
    auto thits = theTransformer->getTransientRecHits(ttrack);

    std::cout << "track resulted in " << trajs.size() << " trajectories and " << thits.size() << " hits!" << std::endl;
    for( const auto& trj : trajs ) {
      std::cout << "original track chi2: " << trj.chiSquared() 
		<< " ndof: " << trj.ndof() << std::endl;
    }

    TransientTrackingRecHit::ConstRecHitContainer mtdthits;
    for( auto& ahit : tryBTLLayers(track,hits,geo.product(),magfield.product()) ) {
      mtdthits.push_back(ahit);
    }
    // in the future this should include an intermediate refit before propagating to the ETL
    // for now it is ok
    for( auto& ahit : tryETLLayers(track,hits,geo.product(),magfield.product()) ) {
      mtdthits.push_back(ahit);
    }

    std::cout << "got " << mtdthits.size() << " new transient hits" << std::endl;
    
    auto ordering = checkRecHitsOrdering(thits);
    if( ordering == RefitDirection::insideOut) {
      std::cout << "fit is inside-out" << std::endl;
      for( auto& ahit : mtdthits ) thits.push_back(ahit);    
    } else {
      std::cout << "fit is outside-in" << std::endl;
      std::reverse(mtdthits.begin(),mtdthits.end());
      for( auto& ahit : thits ) mtdthits.push_back(ahit);
      thits.swap(mtdthits);
    }
    std::cout << "refitting the track with " << thits.size() << " MTD rechits included" << std::endl;
    auto trajwithmtd = theTransformer->transform(ttrack,thits);
    std::cout << "refitting resulted in " << trajwithmtd.size() << " trajectories!" << std::endl;
    for( const auto& trj : trajwithmtd ) {
      std::cout << "mtd track chi2: " << trj.chiSquared() 
		<< " ndof: " << trj.ndof() << std::endl;
      const auto& thetrj = (updateTraj_ ? trj : trajs.front());
      reco::Track result = buildTrack(track, thetrj, bs, magfield.product(), !mtdthits.empty());
      if( result.ndof() > 0 ) {
	/// setup the track extras
	reco::TrackExtra::TrajParams trajParams;
	reco::TrackExtra::Chi2sFive chi2s; 
	size_t hitsstart = outhits->size();
	if( updatePattern_ ) {
	  t2t(trj,*outhits,trajParams,chi2s); // this fills the output hit collection
	} else {
	  t2t(thetrj,*outhits,trajParams,chi2s);
	}
	size_t hitsend = outhits->size();
	extras->push_back(buildTrackExtra(trj : thetrj));
	extras->back().setHits(hitsRefProd,hitsstart,hitsend);	
	extras->back().setTrajParams(trajParams,chi2s);
	//create the track
	output->push_back(result);
	auto& backtrack = output->back();
	reco::TrackExtraRef extraRef(extrasRefProd,extras->size()-1);
	backtrack.setExtra( (updateExtra_ ? extraRef : track.extra()) );
	for(unsigned ihit = hitsstart; ihit < hitsend; ++ihit) {
	  backtrack.appendHitPattern((*outhits)[ihit],ttopo);
	}
	std::cout << "input track: " 
		  << track.pt() << ' ' << track.eta() << ' ' << track.phi() << ' '
		  << track.chi2() << ' ' << track.ndof() << ' ' << track.dxy() << ' '
		  << track.dz() << std::endl;
	std::cout << "added track: " 
		  << result.pt() << ' ' << result.eta() << ' ' << result.phi() << ' '
		  << result.chi2() << ' ' << result.ndof() << ' ' << result.dxy() << ' '
		  << result.dz() << std::endl;
      }
    }
  }

  std::cout << "from " << tracks.size() << " input tracks --->>> got " << output->size() << " tracks with MTD!" << std::endl;

  ev.put(std::move(output));
  ev.put(std::move(extras));
  ev.put(std::move(outhits));
}

template<class TrackCollection>
TransientTrackingRecHit::ConstRecHitContainer
TrackExtenderWithMTDT<TrackCollection>::tryBTLLayers(const TrackType& track,
						     const MTDTrackingDetSetVector& hits,
						     const MTDDetLayerGeometry* geo,
						     const MagneticField* field) const {
  TransientTrackingRecHit::ConstRecHitContainer output;
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
    
    PropagatorWithMaterial prop(anyDirection,0.13957018,field,1.6,false,0.1,true);

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
	  LocalError err = detWithState.second.localError().positionError();
	  
	  cout << "compatibleDets: " << compDets.size() << endl
	       << "  final state pos: " << prop_pos << endl 
	       << "  local state pos: " << local_pos << ' ' << std::sqrt(err.xx()) << ' ' << std::sqrt(err.yy()) << endl
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
	    auto best = detitr->end();
	    double best_chi2 = std::numeric_limits<double>::max();
	    for( auto itr = detitr->begin(); itr != detitr->end(); ++itr ) {
	      auto est =  theEstimator->estimate(detWithState.second,*itr);
	      std::cout << itr->localPosition() << ' ' << itr->localPositionError() << ' ' 
			<< est.first << ' ' << est.second << ' ' << track.chi2() << ' ' << track.ndof() << std::endl;
	      if( est.first && est.second < best_chi2 ) { // just take the best chi2
		best = itr;
		best_chi2 = est.second;
	      }
	    }
	    if( best != detitr->end() ) {
	      output.push_back(hitbuilder->build(&*best));
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
TransientTrackingRecHit::ConstRecHitContainer
TrackExtenderWithMTDT<TrackCollection>::tryETLLayers(const TrackType& track, 
						     const MTDTrackingDetSetVector& hits,
						     const MTDDetLayerGeometry* geo,
						     const MagneticField* field) const {
  TransientTrackingRecHit::ConstRecHitContainer output;
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

    PropagatorWithMaterial prop(anyDirection,0.13957018,field,1.6,false,0.1,true);

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
	  LocalError err = detWithState.second.localError().positionError();
	  
	  cout << "compatibleDets: " << compDets.size() << endl
	       << "  final state pos: " << prop_pos << endl 
	       << "  local state pos: " <<  local_pos << ' ' << std::sqrt(err.xx()) << ' ' << std::sqrt(err.yy()) << endl
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
	    auto best = detitr->end();
	    double best_chi2 = std::numeric_limits<double>::max();
	    for( auto itr = detitr->begin(); itr != detitr->end(); ++itr ) {
	      auto est =  theEstimator->estimate(detWithState.second,*itr);
	      std::cout << itr->localPosition() << ' ' << itr->localPositionError() << ' ' 
			<< est.first << ' ' << est.second << ' ' << track.chi2() << ' ' << track.ndof() << std::endl;
	      if( est.first && est.second < best_chi2 ) { // just take the best chi2
		best = itr;
		best_chi2 = est.second;
	      }
	    }
	    if( best != detitr->end() ) {
	      output.push_back(hitbuilder->build(&*best));
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


//below is unfortunately ripped from other places but 
//since track producer doesn't know about MTD we have to do this
template<class TrackCollection>
reco::Track TrackExtenderWithMTDT<TrackCollection>::buildTrack(const reco::Track& orig,
							       const Trajectory& traj,
							       const reco::BeamSpot& bs,
							       const MagneticField* field,
							       bool hasMTD) const {  
  // get the state closest to the beamline
  TrajectoryStateOnSurface stateForProjectionToBeamLineOnSurface = 
    traj.closestMeasurement(GlobalPoint(bs.x0(),bs.y0(),bs.z0())).updatedState();
  
  if UNLIKELY(!stateForProjectionToBeamLineOnSurface.isValid()){
    edm::LogError("CannotPropagateToBeamLine")<<"the state on the closest measurement isnot valid. skipping track.";
    return reco::Track();
  }

  const FreeTrajectoryState & stateForProjectionToBeamLine=*stateForProjectionToBeamLineOnSurface.freeState();

  PropagatorWithMaterial thePropagator(anyDirection,0.13957018,field,1.6,false,0.1,true);
  
  TSCBLBuilderWithPropagator tscblBuilder(thePropagator);
  TrajectoryStateClosestToBeamLine tscbl = tscblBuilder(stateForProjectionToBeamLine,bs);
  
  if UNLIKELY(!tscbl.isValid()) {
    return reco::Track();
  }

  GlobalPoint v = tscbl.trackStateAtPCA().position();
  math::XYZPoint  pos( v.x(), v.y(), v.z() );
  GlobalVector p = tscbl.trackStateAtPCA().momentum();
  math::XYZVector mom( p.x(), p.y(), p.z() );

  int ndof = traj.ndof();
  
  return reco::Track(traj.chiSquared(),
		     int(ndof),//FIXME fix weight() in TrackingRecHit
		     pos, mom, tscbl.trackStateAtPCA().charge(), 
		     tscbl.trackStateAtPCA().curvilinearError(),
		     orig.algo(),reco::TrackBase::undefQuality,0,0,hasMTD ? 1.0 : -1.0,hasMTD ? 1.0 : -1.0);
}

template<class TrackCollection>
reco::TrackExtra TrackExtenderWithMTDT<TrackCollection>::buildTrackExtra(const Trajectory& trajectory) const {

  const string metname = "MTD|RecoMTD|TrackExtenderWithMTD";

  const Trajectory::RecHitContainer transRecHits = trajectory.recHits();
  
  // put the collection of TrackingRecHit in the event
  
  // sets the outermost and innermost TSOSs
  // FIXME: check it!
  TrajectoryStateOnSurface outerTSOS;
  TrajectoryStateOnSurface innerTSOS;
  unsigned int innerId=0, outerId=0;
  TrajectoryMeasurement::ConstRecHitPointer outerRecHit;
  DetId outerDetId;

  if (trajectory.direction() == alongMomentum) {
    LogTrace(metname)<<"alongMomentum";
    outerTSOS = trajectory.lastMeasurement().updatedState();
    innerTSOS = trajectory.firstMeasurement().updatedState();
    outerId = trajectory.lastMeasurement().recHit()->geographicalId().rawId();
    innerId = trajectory.firstMeasurement().recHit()->geographicalId().rawId();
    outerRecHit =  trajectory.lastMeasurement().recHit();
    outerDetId =   trajectory.lastMeasurement().recHit()->geographicalId();
  } 
  else if (trajectory.direction() == oppositeToMomentum) {
    LogTrace(metname)<<"oppositeToMomentum";
    outerTSOS = trajectory.firstMeasurement().updatedState();
    innerTSOS = trajectory.lastMeasurement().updatedState();
    outerId = trajectory.firstMeasurement().recHit()->geographicalId().rawId();
    innerId = trajectory.lastMeasurement().recHit()->geographicalId().rawId();
    outerRecHit =  trajectory.firstMeasurement().recHit();
    outerDetId =   trajectory.firstMeasurement().recHit()->geographicalId();
  }
  else LogError(metname)<<"Wrong propagation direction!";
  
  const GeomDet *outerDet = gtg->idToDet(outerDetId);
  GlobalPoint outerTSOSPos = outerTSOS.globalParameters().position();
  bool inside = outerDet->surface().bounds().inside(outerDet->toLocal(outerTSOSPos));

  
  GlobalPoint hitPos = (outerRecHit->isValid()) ? outerRecHit->globalPosition() :  outerTSOS.globalParameters().position() ;
  
  if(!inside) {
    LogTrace(metname)<<"The Global Muon outerMostMeasurementState is not compatible with the recHit detector! Setting outerMost postition to recHit position if recHit isValid: " << outerRecHit->isValid();
    LogTrace(metname)<<"From " << outerTSOSPos << " to " <<  hitPos;
  }
  
  
  //build the TrackExtra
  GlobalPoint v = (inside) ? outerTSOSPos : hitPos ;
  GlobalVector p = outerTSOS.globalParameters().momentum();
  math::XYZPoint  outpos( v.x(), v.y(), v.z() );   
  math::XYZVector outmom( p.x(), p.y(), p.z() );
  
  v = innerTSOS.globalParameters().position();
  p = innerTSOS.globalParameters().momentum();
  math::XYZPoint  inpos( v.x(), v.y(), v.z() );   
  math::XYZVector inmom( p.x(), p.y(), p.z() );

  reco::TrackExtra trackExtra(outpos, outmom, true, inpos, inmom, true,
                              outerTSOS.curvilinearError(), outerId,
                              innerTSOS.curvilinearError(), innerId,
			      trajectory.direction(),trajectory.seedRef());
  
  return trackExtra;
 
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
//typedef TrackExtenderWithMTDT<reco::GsfTrackCollection> GSFTrackExtenderWithMTD;

DEFINE_FWK_MODULE(TrackExtenderWithMTD);
//DEFINE_FWK_MODULE(GSFTrackExtenderWithMTD);
