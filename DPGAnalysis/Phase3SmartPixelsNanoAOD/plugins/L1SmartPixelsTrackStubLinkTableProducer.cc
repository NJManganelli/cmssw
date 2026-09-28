// Track -> OT stub LINK table: one row per (track, stub) pair across all tracks of a
// TTTrack collection (OT tracklet tracks or SmartPixels tracks), walking getStubRefs():
//   trackIdx  -> row of the track in trackTableName
//   otStubIdx -> row of the stub in the complete OT stub table (L1TOTStub), i.e. its
//                TTStubRef key; the ref must point into `stubs`, the collection L1TOTStub is
//                built from (asserted), so the link is exact and never a float join
// The stub itself (position, layer, bend, truth) is read from L1TOTStub; this table only
// carries the association.

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include <cstdint>
#include <string>
#include <vector>

class L1SmartPixelsTrackStubLinkTableProducer : public edm::global::EDProducer<> {
public:
  using TTTrackCollection = std::vector<TTTrack<Ref_Phase2TrackerDigi_>>;

  explicit L1SmartPixelsTrackStubLinkTableProducer(const edm::ParameterSet& cfg)
      : tracksToken_(consumes<TTTrackCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
        stubsToken_(consumes<TTStubDetSetVec>(cfg.getParameter<edm::InputTag>("stubs"))),
        linkTableName_(cfg.getParameter<std::string>("linkTableName")),
        trackTableName_(cfg.getParameter<std::string>("trackTableName")),
        stubTableName_(cfg.getParameter<std::string>("stubTableName")) {
    produces<nanoaod::FlatTable>();
  }

  void produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup&) const override {
    const auto& tracks = iEvent.get(tracksToken_);
    const auto stubHandle = iEvent.getHandle(stubsToken_);

    std::vector<int32_t> trackIdx, otStubIdx;
    for (size_t it = 0; it < tracks.size(); ++it) {
      for (const auto& stubRef : tracks[it].getStubRefs()) {
        if (stubRef.id() != stubHandle.id() || stubRef.key() >= stubHandle->dataSize())
          throw cms::Exception("SmartPixelsStubProductMismatch")
              << linkTableName_ << ": a track stub ref points to product " << stubRef.id() << " key "
              << stubRef.key() << ", not into the OT stub collection " << stubHandle.id() << " (size "
              << stubHandle->dataSize() << ") that " << stubTableName_ << " is built from";
        trackIdx.push_back(static_cast<int32_t>(it));
        otStubIdx.push_back(static_cast<int32_t>(stubRef.key()));
      }
    }

    auto table = std::make_unique<nanoaod::FlatTable>(trackIdx.size(), linkTableName_, false, false);
    table->addColumn<int32_t>("trackIdx", trackIdx, "row of the track in " + trackTableName_);
    table->addColumn<int32_t>("otStubIdx", otStubIdx, "row of the stub in " + stubTableName_ + " (row == stub key)");
    table->setDoc("OT stubs of the " + trackTableName_ + " tracks (one row per track-stub pair; the stub is " +
                  stubTableName_ + "[otStubIdx])");
    iEvent.put(std::move(table));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("tracks", edm::InputTag("l1tTTTracksFromTrackletEmulation", "Level1TTTracks"));
    desc.add<edm::InputTag>("stubs", edm::InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"))
        ->setComment("the stub collection the stub table is built from; otStubIdx indexes it");
    desc.add<std::string>("linkTableName", "L1TTrackStub");
    desc.add<std::string>("trackTableName", "L1TTrack");
    desc.add<std::string>("stubTableName", "L1TOTStub");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const edm::EDGetTokenT<TTTrackCollection> tracksToken_;
  const edm::EDGetTokenT<TTStubDetSetVec> stubsToken_;
  const std::string linkTableName_;
  const std::string trackTableName_;
  const std::string stubTableName_;
};

DEFINE_FWK_MODULE(L1SmartPixelsTrackStubLinkTableProducer);
