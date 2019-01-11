#include "RecoLocalFastTime/FTLCommonAlgos/interface/MTDTimeCalib.h"

#include "DataFormats/ForwardDetId/interface/BTLDetId.h"
#include "DataFormats/ForwardDetId/interface/ETLDetId.h"

#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/MTDGeometryBuilder/interface/ProxyMTDTopology.h"
#include "Geometry/MTDGeometryBuilder/interface/RectangularMTDTopology.h"
#include "Geometry/MTDNumberingBuilder/interface/GeometricTimingDet.h"

MTDTimeCalib::MTDTimeCalib(edm::ParameterSet const& conf, 
			   const MTDGeometry* geom, const MTDTopology* topo):
  geom_(geom),
  topo_(topo),
  btlTimeOffset_( conf.getParameter<double>("BTLTimeOffset") ),
  etlTimeOffset_( conf.getParameter<double>("ETLTimeOffset") )
{
}
 

float MTDTimeCalib::getTimeCalib(const MTDDetId& id) const
{
  if (id.subDetector() != MTDDetId::FastTime)
    {
      throw cms::Exception("MTDTimeCalib") << "MTDDetId: " << std::hex
						      << id.rawId()
						      << " is invalid!" << std::dec
						      << std::endl;
    }

  float time_calib = 0.;

  if ( id.mtdSubDetector() == MTDDetId::BTL )
    {
      time_calib += btlTimeOffset_;
      BTLDetId hitId(id);
      DetId geoId = hitId.geographicalId( (BTLDetId::CrysLayout) topo_->getMTDTopologyMode() ); //for BTL topology gives different layout id
      const MTDGeomDet* thedet = geom_->idToDet(geoId);
      
      if( thedet == nullptr ) {
	throw cms::Exception("BTLBarDeviceSim") << "GeographicalID: " << std::hex
					      << geoId.rawId()
						<< " (" << id.rawId()<< ") is invalid!" << std::dec
						<< std::endl;
      }
      const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(thedet->topology());
      const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());    
      
      if ( topo_->getMTDTopologyMode() == (int ) BTLDetId::CrysLayout::tile )
	{
	  constexpr float lightCollTime = 0.2;
	  time_calib -= lightCollTime; //simply remove the offset introduced at sim level
	}
      else if ( topo_->getMTDTopologyMode() == (int ) BTLDetId::CrysLayout::bar )
	{
	  //for bar staggered as used in 
	  constexpr float lightSlopeColl = 0.075;
	  time_calib -= 0.5*topo.pitch().first*lightSlopeColl; //time offset for bar time is L/2v 
	}
      else 
	{
	  constexpr float lightSlopeColl = 0.075;
	  time_calib -= 0.5*topo.pitch().second*lightSlopeColl; //time offset for bar time is L/2v 
	}
    }
  else if ( id.mtdSubDetector() == MTDDetId::ETL )
    {
      time_calib += etlTimeOffset_;
    }
  else
    {
      throw cms::Exception("MTDThresholdClusterizer") << "MTDDetId: " << std::hex
						      << id.rawId()
						      << " is invalid!" << std::dec
						      << std::endl;
    }

  return time_calib;
}

#include "FWCore/Utilities/interface/typelookup.h"

//--- Now use the Framework macros to set it all up:
TYPELOOKUP_DATA_REG(MTDTimeCalib);
