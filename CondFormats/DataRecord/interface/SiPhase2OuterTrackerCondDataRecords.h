#ifndef CondFormats_SiPhase2OuterTrackerCondDataRecords_h
#define CondFormats_SiPhase2OuterTrackerCondDataRecords_h

#include "FWCore/Framework/interface/EventSetupRecordImplementation.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "FWCore/Utilities/interface/mplVector.h"

/*Record associated to SiPhase2OuterTrackerLorentzAngle Object: the SimRcd is used in simulation only*/
class SiPhase2OuterTrackerLorentzAngleRcd
    : public edm::eventsetup::DependentRecordImplementation<SiPhase2OuterTrackerLorentzAngleRcd,
                                                            edm::mpl::Vector<IdealGeometryRecord, TrackerTopologyRcd> > {
};
class SiPhase2OuterTrackerLorentzAngleSimRcd
    : public edm::eventsetup::DependentRecordImplementation<SiPhase2OuterTrackerLorentzAngleSimRcd,
                                                            edm::mpl::Vector<IdealGeometryRecord, TrackerTopologyRcd> > {
};
/*Record associated to SiStripBadStrip Object:*/
class SiPhase2OuterTrackerBadStripRcd : public edm::eventsetup::DependentRecordImplementation<
                                            SiPhase2OuterTrackerBadStripRcd,
                                            edm::mpl::Vector<TrackerTopologyRcd, TrackerDigiGeometryRecord> > {};
#endif
