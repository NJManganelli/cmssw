# Procurement keep-list for the next 200PU SmartPixels inputs

Purpose: the precise per-event content the new 200PU RelVal/production inputs must
persist so that BOTH of the following run in-job with NO product-not-found:

1. **Posture C** — rebuild new-layout L1 tracks from the file's stub tier (already
   documented in `PostureGapStudy.md`; reproduced compactly in §1 here).
2. **PF / Puppi / jet re-emulation** — refill the SC4/SC8 jet + jet-constituent
   link nano tables (`L1puppiJetSC{4,8}`, `L1SC{4,8}JetCands`, `L1PuppiCand`) in-job.

This feeds the procurement request directly. It is derived from the reduced-menu
inventory of the current 100-evt PU RelVal
(`RelValTTbar_14TeV_PU_150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU-v1.root`,
inventory `spxsmoke/pu_file_labels.txt`), which is MISSING the deregionized Puppi
collection and forced the SmartPix nano to drop the PuppiCand/jet-link family. See
§3 for what the current file does vs. does not permit.

---

## 1. Posture-C core (track side) — unchanged from PostureGapStudy.md

Required products (process label shown as it appears in the current PU file, HLT):

- `TTStubsFromPhase2TrackerDigis` : `StubAccepted` (+ `ClusterAccepted`); the tracklet
  chain consumes STUBS, so tracks are rebuilt without re-DIGI.
- Truth maps (Accepted set at minimum):
  `TTClusterAssociatorFromPixelDigis`, `TTStubAssociatorFromPixelDigis`,
  `TTTrackAssociatorFromPixelDigis` (+`…Extended` if extended tracks are wanted).
- `mix:MergedTrackTruth` TrackingParticles **and** TrackingVertices
  (TP::vz() dereferences the parent TrackingVertex — keep both).
- `offlineBeamSpot` (tracklet fit + associator input).

digiRefit (IT) inputs, additionally:
- `simSiPixelDigis:Pixel` digis + PixelDigiSimLinks.
- `g4SimHits:TrackerHitsPixelBarrelLowTof` (PSimHit projector diagnostic) and
  `g4SimHits` SimTracks (parent-momentum angle synthesis).

With these, the in-job chain
`ProducerDTC → l1tTTTracksFromTrackletEmulation → TTTrackAssociatorFromPixelDigis`
yields new-layout tracks with real helix covariance and a fresh truth map
(`seedCovMode=trackCov` valid on PU).

---

## 2. PF / Puppi / jet re-emulation (jet side) — NEW requirement this doc adds

The SC4/SC8 jets and their constituent-link tables need the **flat deregionized
Puppi candidate collection** and the jets **clustered from that same collection
in-job**. Two acceptable procurement options; pick ONE per collection family.

### Option A (preferred, smallest on disk): persist the deregionizer INPUT, re-emulate in-job
Keep only the per-region Puppi PF candidates; re-run the deregionizer + SeededCone
jets in the SmartPix job (helper `reemulateJetSideForPFTier` in
`DPGAnalysis/Phase3SmartPixelsNanoAOD/python/l1tPh3SmartPixelsNano_cff.py`):

REQUIRED product:
- `l1tLayer1:PuppiRegional`  (`l1t::RegionalOutput<vector<l1t::PFCandidate>>`)
  — the DeregionizerProducer input.  (Barrel/HGCal/HF regions all come from this
  single product in the current menu.)

In-job chain (all downstream products regenerated, so their Ptrs are self-consistent):
```
l1tLayer1:PuppiRegional  →  DeregionizerProducer          (l1tLayer2Deregionizer:Puppi)
                         →  L1SeedConePFJetProducer SC4   (l1tSC4PFL1PuppiCorrectedEmulator)
                         →  L1SeedConePFJetProducer SC8   (l1tSC8PFL1PuppiCorrectedEmulator, coneSize 0.8)
```
`L1PuppiCand.l1TrackIdx` crossref additionally needs the posture-C
`l1tTTTracksFromTrackletEmulation` (already in §1).

This is the current file's ONLY viable route IF `l1tLayer1:PuppiRegional` is kept
(it is present in the current file — see §3).

### Option B (fallback, larger): persist the deregionizer OUTPUT and the jets together
Keep both, as one coherent set (jets Ptr into the kept Puppi collection):
- `l1tLayer2Deregionizer:Puppi`            (`vector<l1t::PFCandidate>`)
- `l1tSC4PFL1PuppiCorrectedEmulator`       (`vector<l1t::PFJet>`)
- `l1tSC8PFL1PuppiCorrectedEmulator`       (`vector<l1t::PFJet>`)
- (optional, extended/displaced tier) `l1tLayer2DeregionizerExtended:Puppi`,
  `l1tSC4PFL1PuppiExtendedCorrectedEmulator`.

**Coherence rule (the trap the current file falls into):** if the jets are kept but
the deregionized Puppi they were clustered from is NOT (or is a different process),
every `candIdx` in the link table resolves to `-1` (dangling Ptr). Keep the jets
and their source Puppi from the SAME process, or keep only the regional input and
re-emulate (Option A). Do not mix.

### Optional jet-side extras (only if those tables are wanted)
- NG SC4 jets: `l1tSC4NGJetProducer:l1tSC4NGJets` + `l1tLayer2DeregionizerExtended:Puppi`
  + `l1tTTTracksFromExtendedTrackletEmulation` (for `L1SC4NGJetCands`/`L1ExtPuppiCand`).
  NOT needed for the SC8 training ask (SC8 uses the PLAIN deregionizer Puppi + plain
  tracks; there is no NG SC8 producer until a model is deployed).
- HGCal cluster crossref (`L1HGCCluster`, `hgcClusterIdx`):
  `l1tHGCalBackEndLayer2Producer:HGCalBackendLayer2Processor3DClustering`
  (present in the current file).

---

## 3. What the CURRENT 100-evt PU file permits (feasibility verdict)

From `spxsmoke/pu_file_labels.txt` + `edmDumpEventContent`:

| Need                                   | In current file?                          |
|----------------------------------------|-------------------------------------------|
| `l1tLayer1:PuppiRegional`              | **YES** (regional Puppi PFCandidates)     |
| `l1tLayer2Deregionizer:Puppi` (flat)   | **NO** (the reduced-menu gap)             |
| `l1tSC8PFL1PuppiCorrectedEmulator`     | YES — but constituents dangle (no source) |
| `l1tSC4PFL1PuppiCorrectedEmulator`     | YES — same dangling caveat                |
| posture-C stub tier + maps + TPs + BS  | YES (all present)                         |
| `l1tHGCalBackEndLayer2Producer` clusters | YES                                      |

**Verdict: the jet-side chain IS re-runnable in-job on the CURRENT file via Option A**
(re-emulate the deregionizer + SC4/SC8 jets from the persisted `l1tLayer1:PuppiRegional`).
Consuming the file's pre-persisted SC8 jets directly does NOT work — their
constituents Ptr into the absent file deregionizer product, so `candIdx` would be
all `-1`. The SmartPix nano therefore should re-emulate rather than consume file jets.

Because the current file already carries `l1tLayer1:PuppiRegional`, jet-side studies
can proceed on it now. The new procurement only needs to GUARANTEE that this product
(Option A) — or the coherent Option-B set — is in the kept event content, so the
gap does not recur.

---

## 4. Minimal procurement keep block (Option A, recommended)

```
# --- posture-C track side ---
keep *_TTStubsFromPhase2TrackerDigis_*_*
keep *_TTClusterAssociatorFromPixelDigis_*_*
keep *_TTStubAssociatorFromPixelDigis_*_*
keep *_TTTrackAssociatorFromPixelDigis_*_*
keep TrackingParticles_mix_MergedTrackTruth_*
keep TrackingVertexs_mix_MergedTrackTruth_*
keep *_offlineBeamSpot_*_*
# --- digiRefit IT side ---
keep PixelDigiSimLink*_simSiPixelDigis_Pixel_*
keep *_simSiPixelDigis_Pixel_*
keep *_g4SimHits_TrackerHitsPixelBarrelLowTof_*
keep SimTracks_g4SimHits_*_*
# --- PF/Puppi/jet side (Option A: re-emulate from regional Puppi) ---
keep *_l1tLayer1_PuppiRegional_*
# optional endcap cluster crossref
keep *_l1tHGCalBackEndLayer2Producer_HGCalBackendLayer2Processor3DClustering_*
```

If Option B is chosen instead of Option A, replace the single
`l1tLayer1_PuppiRegional` keep with the four `l1tLayer2Deregionizer:Puppi` +
`l1tSC{4,8}PFL1PuppiCorrectedEmulator` keeps from §2 Option B (jets and their source
Puppi from the same process).
