// Emits an extension table on a L1 track table carrying the track fit's helix
// covariance, TTTrack::helixCovMat(): the 15 independent elements of the
// symmetric 5x5 matrix, upper triangle, named cov_<p>_<q>.
//
// PARAMETER ORDER AND UNITS are TTTrack's own (enum Hpar): rInv [1/cm],
// phi [rad], tanL [unitless], z0 [cm], d0 [cm] -- the same quantities as the
// main table's rInv/phi/tanL/z0/d0 columns, in the same sign conventions, so
// each column's square root is the fit's claimed uncertainty on that column.
// For the hybrid chain the matrix is the TMTT Kalman state covariance with the
// 1/2R row and column rescaled to 1/R (KFParamsComb::trackParamsCov). The
// azimuth in the KF state is sector-local; a constant rotation leaves the
// covariance unchanged.
//
// ALL ZERO when the producer of the tracks did not fill the matrix (a
// default-constructed TTTrack::CovMat): that is "not provided", not "exact".
//
// Index/row alignment relies on the track table using the same source
// collection with an empty cut and native ordering.
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include <array>
#include <string>
#include <vector>

class L1TrackHelixCovTableProducer : public edm::global::EDProducer<> {
public:
  using L1Track = TTTrack<Ref_Phase2TrackerDigi_>;
  using L1TrackCollection = std::vector<L1Track>;

  explicit L1TrackHelixCovTableProducer(const edm::ParameterSet& cfg)
      : tracksToken_(consumes<L1TrackCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
        trackTableName_(cfg.getParameter<std::string>("trackTableName")) {
    produces<nanoaod::FlatTable>();
  }

  void produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup&) const override {
    const auto& tracks = iEvent.get(tracksToken_);
    const unsigned int n = tracks.size();

    // names and units indexed by TTTrack's Hpar enum
    static constexpr std::array<const char*, 5> kName = {{"rInv", "phi", "tanL", "z0", "d0"}};
    static constexpr std::array<const char*, 5> kUnit = {{"1/cm", "rad", "1", "cm", "cm"}};
    static_assert(L1Track::INVR == 0 && L1Track::PHI0 == 1 && L1Track::TANL == 2 && L1Track::Z0 == 3 &&
                      L1Track::D0 == 4,
                  "column naming assumes TTTrack::Hpar = {INVR, PHI0, TANL, Z0, D0}");

    std::array<std::vector<float>, 15> cols;
    for (auto& c : cols)
      c.resize(n);
    for (unsigned int t = 0; t < n; ++t) {
      const auto& C = tracks[t].helixCovMat();
      unsigned int k = 0;
      for (unsigned int i = 0; i < 5; ++i)
        for (unsigned int j = i; j < 5; ++j)
          cols[k++][t] = C(i, j);
    }

    auto table = std::make_unique<nanoaod::FlatTable>(n, trackTableName_, false, true);
    unsigned int k = 0;
    for (unsigned int i = 0; i < 5; ++i) {
      for (unsigned int j = i; j < 5; ++j) {
        const std::string name = std::string("cov_") + kName[i] + "_" + kName[j];
        const std::string unit = (i == j) ? std::string(kUnit[i]) + "^2"
                                          : std::string(kUnit[i]) + "*" + std::string(kUnit[j]);
        const std::string doc = (i == j) ? std::string("track-fit variance of ") + kName[i] + " [" + unit +
                                               "] (TTTrack::helixCovMat; all-zero covariance = not provided)"
                                         : std::string("track-fit covariance of ") + kName[i] + " and " + kName[j] +
                                               " [" + unit + "] (TTTrack::helixCovMat)";
        // full float precision: off-diagonal terms are small differences of large
        // correlated quantities, and a trimmed mantissa can make the matrix non-PSD
        table->addColumn<float>(name, cols[k++], doc);
      }
    }
    iEvent.put(std::move(table));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("tracks", edm::InputTag("l1tTTTracksFromTrackletEmulation", "Level1TTTracks"));
    desc.add<std::string>("trackTableName", "L1TTrack");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const edm::EDGetTokenT<L1TrackCollection> tracksToken_;
  const std::string trackTableName_;
};

DEFINE_FWK_MODULE(L1TrackHelixCovTableProducer);
