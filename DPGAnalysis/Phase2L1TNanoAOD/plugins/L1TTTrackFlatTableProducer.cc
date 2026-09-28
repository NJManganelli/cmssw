// Flat table of an L1 TTTrack collection: the configured `variables` plus, unless
// helixCov is False, the track fit's helix covariance TTTrack::helixCovMat() as
// the 15 independent elements of the symmetric 5x5 matrix, upper triangle, named
// cov_<p>_<q>. Being part of the track table itself, the covariance always
// describes the table's own rows (same src, same cut).
//
// PARAMETER ORDER AND UNITS are TTTrack's own (enum Hpar): rInv [1/cm],
// phi [rad], tanL [unitless], z0 [cm], d0 [cm] -- the same quantities as the
// rInv/phi/tanL/z0/d0 columns, in the same sign conventions, so each column's
// square root is the fit's claimed uncertainty on that column. For the hybrid
// chain the matrix is the TMTT Kalman state covariance with the 1/2R row and
// column rescaled to 1/R (KFParamsComb::trackParamsCov). The azimuth in the KF
// state is sector-local; a constant rotation leaves the covariance unchanged.
//
// ALL ZERO when the producer of the tracks did not fill the matrix (a
// default-constructed TTTrack::CovMat): that is "not provided", not "exact".
#include "PhysicsTools/NanoAOD/interface/SimpleFlatTableProducer.h"

#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"

#include <array>
#include <string>
#include <vector>

class L1TTTrackFlatTableProducer : public SimpleFlatTableProducer<TTTrack<Ref_Phase2TrackerDigi_>> {
public:
  using L1Track = TTTrack<Ref_Phase2TrackerDigi_>;
  using Base = SimpleFlatTableProducer<L1Track>;

  explicit L1TTTrackFlatTableProducer(const edm::ParameterSet& params)
      : Base(params), helixCov_(params.getParameter<bool>("helixCov")) {}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc = Base::baseDescriptions();
    desc.add<bool>("helixCov", true)->setComment("add the cov_<p>_<q> helix covariance columns");
    descriptions.addWithDefaultLabel(desc);
  }

  std::unique_ptr<nanoaod::FlatTable> fillTable(const edm::Event& iEvent,
                                                const edm::Handle<edm::View<L1Track>>& prod) const override {
    auto out = Base::fillTable(iEvent, prod);
    if (!helixCov_)
      return out;

    // the rows the base selected, by the same rule
    std::vector<const L1Track*> rows;
    if (prod.isValid() || !this->skipNonExistingSrc_) {
      if (this->singleton_) {
        rows.push_back(&(*prod)[0]);
      } else {
        for (const auto& obj : *prod) {
          if (this->cut_(obj))
            rows.push_back(&obj);
          if (rows.size() >= this->maxLen_)
            break;
        }
      }
    }

    static constexpr std::array<const char*, 5> kName = {{"rInv", "phi", "tanL", "z0", "d0"}};
    static constexpr std::array<const char*, 5> kUnit = {{"1/cm", "rad", "1", "cm", "cm"}};
    static_assert(L1Track::INVR == 0 && L1Track::PHI0 == 1 && L1Track::TANL == 2 && L1Track::Z0 == 3 &&
                      L1Track::D0 == 4,
                  "column naming assumes TTTrack::Hpar = {INVR, PHI0, TANL, Z0, D0}");
    std::vector<float> col(rows.size());
    for (unsigned int i = 0; i < 5; ++i) {
      for (unsigned int j = i; j < 5; ++j) {
        for (size_t t = 0; t < rows.size(); ++t)
          col[t] = rows[t]->helixCovMat()(i, j);
        const std::string name = std::string("cov_") + kName[i] + "_" + kName[j];
        const std::string unit =
            (i == j) ? std::string(kUnit[i]) + "^2" : std::string(kUnit[i]) + "*" + std::string(kUnit[j]);
        const std::string doc = (i == j) ? std::string("track-fit variance of ") + kName[i] + " [" + unit +
                                               "] (TTTrack::helixCovMat; all-zero covariance = not provided)"
                                         : std::string("track-fit covariance of ") + kName[i] + " and " + kName[j] +
                                               " [" + unit + "] (TTTrack::helixCovMat)";
        // full float precision: off-diagonal terms are small differences of large
        // correlated quantities, and a trimmed mantissa can make the matrix non-PSD
        out->addColumn<float>(name, col, doc);
      }
    }
    return out;
  }

private:
  const bool helixCov_;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1TTTrackFlatTableProducer);
