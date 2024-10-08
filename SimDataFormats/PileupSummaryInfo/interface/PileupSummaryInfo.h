#ifndef PileupSummaryInfo_h
#define PileupSummaryInfo_h
// -*- C++ -*-
//
// Package:     PileupSummaryInfo
// Class  :     PileupSummaryInfo
//
/**\class PileupSummaryInfo PileupSummaryInfo.h SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h

Description: contains information related to the details of the pileup simulation for a given event
Usage: purely descriptive
*/
//
// Original Author:  Mike Hildreth, Notre Dame
//         Created:  July 1, 2010
//

#include "DataFormats/Provenance/interface/EventID.h"
#include <vector>
#include <string>

class PileupSummaryInfo {
public:
  PileupSummaryInfo() {}

  PileupSummaryInfo(const int num_PU_vertices,
                    const std::vector<float>& zpositions,
                    const std::vector<float>& times,  // may be empty
                    const std::vector<float>& sumpT_lowpT,
                    const std::vector<float>& sumpT_highpT,
                    const std::vector<int>& ntrks_lowpT,
                    const std::vector<int>& ntrks_highpT);

  PileupSummaryInfo(const int num_PU_vertices,
                    const std::vector<float>& zpositions,
                    const std::vector<float>& times,  // may be empty
                    const std::vector<float>& sumpT_lowpT,
                    const std::vector<float>& sumpT_highpT,
                    const std::vector<int>& ntrks_lowpT,
                    const std::vector<int>& ntrks_highpT,
                    int bunchCrossing);

  PileupSummaryInfo(const int num_PU_vertices,
                    const std::vector<float>& zpositions,
                    const std::vector<float>& times,  // may be empty
                    const std::vector<float>& sumpT_lowpT,
                    const std::vector<float>& sumpT_highpT,
                    const std::vector<int>& ntrks_lowpT,
                    const std::vector<int>& ntrks_highpT,
                    const std::vector<edm::EventID>& eventInfo,
                    const std::vector<float>& pT_hats,
                    int bunchCrossing,
                    float TrueNumInteractions,
                    int bunchSpacing);

  PileupSummaryInfo(const int num_PU_vertices,
                    const std::vector<float>& instLumi,
                    const std::vector<edm::EventID>& eventInfo);

  ~PileupSummaryInfo();

  const int getPU_NumInteractions() const { return num_PU_vertices_; }
  const std::vector<float>& getPU_zpositions() const { return zpositions_; }
  bool has_times() const { return !times_.empty(); }
  const std::vector<float>& getPU_times() const { return times_; }
  const std::vector<float>& getPU_sumpT_lowpT() const { return sumpT_lowpT_; }
  const std::vector<float>& getPU_sumpT_highpT() const { return sumpT_highpT_; }
  const std::vector<int>& getPU_ntrks_lowpT() const { return ntrks_lowpT_; }
  const std::vector<int>& getPU_ntrks_highpT() const { return ntrks_highpT_; }
  const std::vector<float>& getPU_instLumi() const { return instLumi_; }
  const std::vector<edm::EventID>& getPU_EventID() const { return eventInfo_; }
  const std::vector<float>& getPU_pT_hats() const { return pT_hats_; }
  const int getBunchCrossing() const { return bunchCrossing_; }
  const int getBunchSpacing() const { return bunchSpacing_; }
  const float getTrueNumInteractions() const { return TrueNumInteractions_; }

private:
  // for "standard" pileup: we have MC Truth information for these

  int num_PU_vertices_;
  std::vector<float> zpositions_;
  std::vector<float> times_;
  std::vector<float> sumpT_lowpT_;
  std::vector<float> sumpT_highpT_;
  std::vector<int> ntrks_lowpT_;
  std::vector<int> ntrks_highpT_;
  std::vector<edm::EventID> eventInfo_;
  std::vector<float> pT_hats_;
  int bunchCrossing_;
  int bunchSpacing_;
  float TrueNumInteractions_;

  // for DataMixer pileup, we only have raw information:

  std::vector<float> instLumi_;
};

#endif
