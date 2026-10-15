// Nano producer for the DTC stub frames the L1 tracklet track finder reads: one row per (stub, TFP
// region) in ProducerDTC's output streams (TTDTC, ProducerDTC:StubAccepted).
//
// WHY THIS EXISTS. L1TOTStub describes each stub geometrically (float position, FE bend). The
// tracklet chain does not read that: L1FPGATrackProducer takes, per processing region (nonant)
// and DTC channel, the DTC's 64-bit frame, whose 36-bit stub word carries the integer r / z / phi,
// the back-end bend code and the layer code (L1FPGATrackProducer.cc, "READ DTC STUB INFORMATION").
// Emulating the finder from L1TOTStub alone misses four things this table gives exactly:
//   - the integer coordinates and the bend code as the finder decodes them (trklet::Stub);
//   - the nonant routing, including the overlap copy of a stub near a nonant border;
//   - DTC acceptance (a stub absent from every stream was dropped by the DTC, e.g. by truncation);
//   - the stream order, which fixes the order stubs reach the input router and the VM memories.
//
// ROWS. Every non-null frame of every (region, channel) stream, in the loop order of
// L1FPGATrackProducer (region, then channel, then position in the stream). Gaps (null frames) are
// not rows; `index` keeps the frame's position in its stream, gaps included.
//
// otStubIdx is the L1TOTStub row of the frame's TTStubRef (that table is in collection order, row ==
// TTStubRef key). The producer asserts the DTC consumed the same stub product as L1TOTStub.
//
// withPositions (off by default; the closure's one-file checks) adds the TTStub position the
// tracklet KF fits (trackerDTC::Setup::stubPosTT) and the bit-accurate position decoded from the
// frame (Setup::stubPosDTC), with the phi LSB of the frame's module type.

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/L1TrackTrigger/interface/TTDTC.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "L1Trigger/TrackerDTC/interface/Setup.h"
#include "L1Trigger/TrackerDTC/interface/SetupRcd.h"

#include <cstdint>
#include <string>
#include <vector>

class L1SmartPixelsDTCStubTableProducer : public edm::stream::EDProducer<> {
public:
  explicit L1SmartPixelsDTCStubTableProducer(const edm::ParameterSet& cfg)
      : dtcToken_(consumes<TTDTC>(cfg.getParameter<edm::InputTag>("dtc"))),
        stubsToken_(consumes<TTStubDetSetVec>(cfg.getParameter<edm::InputTag>("stubs"))),
        tableName_(cfg.getParameter<std::string>("tableName")),
        withPositions_(cfg.getParameter<bool>("withPositions")) {
    if (withPositions_)
      setupToken_ = esConsumes<trackerDTC::Setup, trackerDTC::SetupRcd>();
    produces<nanoaod::FlatTable>();
  }

  void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override {
    const TTDTC& dtc = iEvent.get(dtcToken_);
    edm::Handle<TTStubDetSetVec> stubHandle;
    iEvent.getByToken(stubsToken_, stubHandle);
    const trackerDTC::Setup* setup = withPositions_ ? &iSetup.getData(setupToken_) : nullptr;

    std::vector<int32_t> otStubIdx;
    std::vector<uint8_t> region, channel;
    std::vector<uint16_t> index;
    std::vector<uint32_t> frameHi, frameLo;
    std::vector<float> xTT, yTT, zTT, rDTC, phiDTC, zDTC, phiBase;

    for (const int reg : dtc.tfpRegions()) {
      for (const int chan : dtc.tfpChannels()) {
        const tt::StreamStub& stream = dtc.stream(reg, chan);
        for (size_t i = 0; i < stream.size(); ++i) {
          const tt::FrameStub& frame = stream[i];
          if (frame.first.isNull())
            continue;
          if (frame.first.id() != stubHandle.id() || frame.first.key() >= stubHandle->dataSize())
            throw cms::Exception("SmartPixelsDTCStubProductMismatch")
                << tableName_ << ": a DTC frame refers to stub product " << frame.first.id() << " key "
                << frame.first.key() << ", not to the L1TOTStub product " << stubHandle.id()
                << "; otStubIdx would point at the wrong stub";
          const unsigned long long word = frame.second.to_ullong();
          otStubIdx.push_back(static_cast<int32_t>(frame.first.key()));
          region.push_back(static_cast<uint8_t>(reg));
          channel.push_back(static_cast<uint8_t>(chan));
          index.push_back(static_cast<uint16_t>(i));
          frameHi.push_back(static_cast<uint32_t>(word >> 32));
          frameLo.push_back(static_cast<uint32_t>(word & 0xffffffffULL));
          if (withPositions_) {
            const GlobalPoint tt = setup->stubPosTT(frame.first);
            const GlobalPoint hw = setup->stubPosDTC(frame, reg);
            xTT.push_back(tt.x());
            yTT.push_back(tt.y());
            zTT.push_back(tt.z());
            rDTC.push_back(hw.perp());
            phiDTC.push_back(hw.phi());
            zDTC.push_back(hw.z());
            phiBase.push_back(setup->stubBasePhi(setup->sensorModule(frame.first)->type()));
          }
        }
      }
    }

    auto tab = std::make_unique<nanoaod::FlatTable>(otStubIdx.size(), tableName_, false, false);
    tab->addColumn<int32_t>("otStubIdx", otStubIdx, "row of the stub in L1TOTStub (row == TTStubRef key)");
    tab->addColumn<uint8_t>("region", region, "TFP processing region (nonant) 0-8 of this stream");
    tab->addColumn<uint8_t>("channel", channel, "TFP input channel 0-47 (DTC link) of this stream");
    tab->addColumn<uint16_t>("index", index, "position of the frame in its stream, gaps included");
    tab->addColumn<uint32_t>("frameHi", frameHi,
                             "DTC frame (tt::Frame, 64 bits) bits 63-32. Layout (StubDTC.cc, L1FPGATrackProducer.cc): "
                             "bit 0 valid, bits 2-1 layer code, bits 38-3 the 36-bit stub word the tracklet chain "
                             "decodes (trklet::Stub: from the top, negDisk, r, z, phi, alpha, bend)");
    tab->addColumn<uint32_t>("frameLo", frameLo, "DTC frame bits 31-0");
    if (withPositions_) {
      tab->addColumn<float>("xTT", xTT, "trackerDTC::Setup::stubPosTT x [cm]: the position the tracklet KF fits");
      tab->addColumn<float>("yTT", yTT, "trackerDTC::Setup::stubPosTT y [cm]");
      tab->addColumn<float>("zTT", zTT, "trackerDTC::Setup::stubPosTT z [cm]");
      tab->addColumn<float>("rDTC", rDTC, "bit-accurate r decoded from the frame (Setup::stubPosDTC) [cm]");
      tab->addColumn<float>("phiDTC", phiDTC, "bit-accurate phi decoded from the frame (Setup::stubPosDTC) [rad]");
      tab->addColumn<float>("zDTC", zDTC, "bit-accurate z decoded from the frame (Setup::stubPosDTC) [cm]");
      tab->addColumn<float>("phiBase", phiBase, "phi LSB of the frame's module type, Setup::stubBasePhi(type) [rad]");
    }
    tab->setDoc(
        "DTC stub frames as the L1 tracklet track finder reads them: one row per (stub, TFP region) of "
        "ProducerDTC:StubAccepted, in the finder's loop order (region, channel, stream position). "
        "otStubIdx is the L1TOTStub row; a stub in no row was not accepted by the DTC.");
    iEvent.put(std::move(tab));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("dtc", edm::InputTag("ProducerDTC", "StubAccepted"))
        ->setComment("the TTDTC the tracklet chain reads (l1tTTTracksFromTrackletEmulation.InputTagTTDTC)");
    desc.add<edm::InputTag>("stubs", edm::InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"))
        ->setComment("the stub product L1TOTStub describes; every frame must refer to it");
    desc.add<std::string>("tableName", "L1TOTStubDTC");
    desc.add<bool>("withPositions", false)
        ->setComment("add xTT/yTT/zTT (stubPosTT), rDTC/phiDTC/zDTC (stubPosDTC) and phiBase per row");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const edm::EDGetTokenT<TTDTC> dtcToken_;
  const edm::EDGetTokenT<TTStubDetSetVec> stubsToken_;
  edm::ESGetToken<trackerDTC::Setup, trackerDTC::SetupRcd> setupToken_;
  const std::string tableName_;
  const bool withPositions_;
};

DEFINE_FWK_MODULE(L1SmartPixelsDTCStubTableProducer);
