#ifndef PF_PU_AssoMapAlgos_h
#define PF_PU_AssoMapAlgos_h

/**\class PF_PU_AssoMap PF_PU_AssoMap.cc CommonTools/RecoUtils/plugins/PF_PU_AssoMap.cc

 Description: Produces a map with association between tracks and their particular most probable vertex with a quality of this association
*/
//
// Original Author:  Matthias Geisler,32 4-B20,+41227676487,
// $Id: PF_PU_AssoMapAlgos.h,v 1.7 2012/11/21 09:45:04 mgeisler Exp $
//
//

#include <string>

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/OneToManyWithQuality.h"
#include "DataFormats/Common/interface/OneToManyWithQualityGeneric.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Math/interface/Point3D.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertex.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertexFwd.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"

//
// constants, enums and typedefs
//

const double kMass = 0.49765;
const double lamMass = 1.11568;

/*limits for the quality criteria*/

const double tw_90 = 1.e-2;
const double tw_70 = 1.e-1;
const double tw_50 = 2.e-1;

const double sec_70 = 5.;
const double sec_50 = 19.;

const double fin_70 = 1.e-1;
const double fin_50 = 3.e-1;

typedef edm::AssociationMap<edm::OneToManyWithQuality<reco::VertexCollection, reco::TrackCollection, int>>
    TrackToVertexAssMap;
typedef edm::AssociationMap<edm::OneToManyWithQuality<reco::TrackCollection, reco::VertexCollection, int>>
    VertexToTrackAssMap;

typedef std::pair<reco::TrackRef, int> TrackQualityPair;
typedef std::vector<TrackQualityPair> TrackQualityPairVector;

typedef std::pair<reco::VertexRef, int> VertexStepPair;

typedef std::pair<reco::VertexRef, TrackQualityPair> VertexTrackQuality;

typedef std::pair<reco::VertexRef, float> VertexPtsumPair;
typedef std::vector<VertexPtsumPair> VertexPtsumVector;

class PF_PU_AssoMapAlgos {
public:
  //dedicated constructor for the algorithms
  PF_PU_AssoMapAlgos(const edm::ParameterSet& iConfig, edm::ConsumesCollector&& iC) : PF_PU_AssoMapAlgos(iConfig, iC) {}
  PF_PU_AssoMapAlgos(const edm::ParameterSet&, edm::ConsumesCollector&);
  // virtual destructor needed when virtual functions declared
  virtual ~PF_PU_AssoMapAlgos() noexcept(false) {}

  //get all needed collections at the beginning
  virtual void GetInputCollections(edm::Event&, const edm::EventSetup&);

  //create the track-to-vertex and vertex-to-track maps in one go
  std::pair<std::unique_ptr<TrackToVertexAssMap>, std::unique_ptr<VertexToTrackAssMap>> createMappings(
      edm::Handle<reco::TrackCollection> trkcollH);

  //create the track to vertex association map
  std::unique_ptr<TrackToVertexAssMap> CreateTrackToVertexMap(edm::Handle<reco::TrackCollection>);

  //create the vertex to track association map
  std::unique_ptr<VertexToTrackAssMap> CreateVertexToTrackMap(edm::Handle<reco::TrackCollection>);

  //function to sort the vertices in the AssociationMap by the sum of (pT - pT_Error)**2
  std::unique_ptr<TrackToVertexAssMap> SortAssociationMap(TrackToVertexAssMap*, edm::Handle<reco::TrackCollection>);

protected:
  //protected functions

  //create helping vertex vector to remove associated vertices
  std::vector<reco::VertexRef> CreateVertexVector(edm::Handle<reco::VertexCollection>);

  //erase one vertex from the vertex vector
  void EraseVertex(std::vector<reco::VertexRef>&, reco::VertexRef);

  //find an association for a certain track
  VertexStepPair FindAssociation(const reco::TrackRef&,
                                 const std::vector<reco::VertexRef>&,
                                 edm::ESHandle<MagneticField>,
                                 edm::ESHandle<GlobalTrackingGeometry>,
                                 edm::Handle<reco::BeamSpot>,
                                 int);

  //get the quality for a certain association
  int DefineQuality(int, int, double);

private:
  // private methods for internal usage

  //function to find the closest vertex in z for a certain track
  static reco::VertexRef FindClosestZ(const reco::TrackRef, const std::vector<reco::VertexRef>&, double tWeight = 0.);

  //function to find the closest vertex in 3D for a certain track
  static reco::VertexRef FindClosest3D(reco::TransientTrack, const std::vector<reco::VertexRef>&, double tWeight = 0.);

  //function to calculate the deltaR between a vector and a vector connecting two points
  static double dR(const math::XYZPoint&, const math::XYZVector&, edm::Handle<reco::BeamSpot>);

  //function to filter the conversion collection
  static std::unique_ptr<reco::ConversionCollection> GetCleanedConversions(edm::Handle<reco::ConversionCollection>,
                                                                           edm::Handle<reco::BeamSpot>,
                                                                           bool);

  //function to find out if the track comes from a gamma conversion
  static bool ComesFromConversion(const reco::TrackRef, const reco::ConversionCollection&, reco::Conversion*);

  static reco::VertexRef FindConversionVertex(const reco::TrackRef,
                                              const reco::Conversion&,
                                              edm::ESHandle<MagneticField>,
                                              edm::ESHandle<GlobalTrackingGeometry>,
                                              edm::Handle<reco::BeamSpot>,
                                              const std::vector<reco::VertexRef>&,
                                              double);

  //function to filter the Kshort collection
  static std::unique_ptr<reco::VertexCompositeCandidateCollection> GetCleanedKshort(
      edm::Handle<reco::VertexCompositeCandidateCollection>, edm::Handle<reco::BeamSpot>, bool);

  //function to filter the Lambda collection
  static std::unique_ptr<reco::VertexCompositeCandidateCollection> GetCleanedLambda(
      edm::Handle<reco::VertexCompositeCandidateCollection>, edm::Handle<reco::BeamSpot>, bool);

  //function to find out if the track comes from a V0 decay
  static bool ComesFromV0Decay(const reco::TrackRef,
                               const reco::VertexCompositeCandidateCollection&,
                               const reco::VertexCompositeCandidateCollection&,
                               reco::VertexCompositeCandidate*);

  static reco::VertexRef FindV0Vertex(const reco::TrackRef,
                                      const reco::VertexCompositeCandidate&,
                                      edm::ESHandle<MagneticField>,
                                      edm::ESHandle<GlobalTrackingGeometry>,
                                      edm::Handle<reco::BeamSpot>,
                                      const std::vector<reco::VertexRef>&,
                                      double);

  //function to filter the nuclear interaction collection
  static std::unique_ptr<reco::PFDisplacedVertexCollection> GetCleanedNI(edm::Handle<reco::PFDisplacedVertexCollection>,
                                                                         edm::Handle<reco::BeamSpot>,
                                                                         bool);

  //function to find out if the track comes from a nuclear interaction
  static bool ComesFromNI(const reco::TrackRef, const reco::PFDisplacedVertexCollection&, reco::PFDisplacedVertex*);

  static reco::VertexRef FindNIVertex(const reco::TrackRef,
                                      const reco::PFDisplacedVertex&,
                                      edm::ESHandle<MagneticField>,
                                      edm::ESHandle<GlobalTrackingGeometry>,
                                      edm::Handle<reco::BeamSpot>,
                                      const std::vector<reco::VertexRef>&,
                                      double);

  //function to find the vertex with the highest TrackWeight for a certain track
  template <typename TREF>
  static reco::VertexRef TrackWeightAssociation(const TREF&, const std::vector<reco::VertexRef>&);

  // ----------member data ---------------------------

  int input_MaxNumAssociations_;

  edm::EDGetTokenT<reco::VertexCollection> token_VertexCollection_;
  edm::Handle<reco::VertexCollection> vtxcollH;

  edm::EDGetTokenT<reco::BeamSpot> token_BeamSpot_;
  edm::Handle<reco::BeamSpot> beamspotH;

  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> token_bField_;
  edm::ESHandle<MagneticField> bFieldH;
  const edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> token_TrackingGeometry_;
  edm::ESHandle<GlobalTrackingGeometry> trackingGeometryH;

  bool input_doReassociation_;
  bool cleanedColls_;

  edm::EDGetTokenT<reco::ConversionCollection> ConversionsCollectionToken_;
  edm::Handle<reco::ConversionCollection> convCollH;
  std::unique_ptr<reco::ConversionCollection> cleanedConvCollP;

  edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> KshortCollectionToken_;
  edm::Handle<reco::VertexCompositeCandidateCollection> vertCompCandCollKshortH;
  std::unique_ptr<reco::VertexCompositeCandidateCollection> cleanedKshortCollP;

  edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> LambdaCollectionToken_;
  edm::Handle<reco::VertexCompositeCandidateCollection> vertCompCandCollLambdaH;
  std::unique_ptr<reco::VertexCompositeCandidateCollection> cleanedLambdaCollP;

  edm::EDGetTokenT<reco::PFDisplacedVertexCollection> NIVertexCollectionToken_;
  edm::Handle<reco::PFDisplacedVertexCollection> displVertexCollH;
  std::unique_ptr<reco::PFDisplacedVertexCollection> cleanedNICollP;

  int input_FinalAssociation_;

  bool ignoremissingpfcollection_;
  bool missingColls;  // is true if there is a diplaced vertex collection in the event

  double input_nTrack_;

  int maxNumWarnings_;  // CV: print Warning if TrackExtra objects don't exist in input file,
  int numWarnings_;     //     but only a few times
};

#endif
