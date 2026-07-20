# L1Trigger/Phase3SmartPixels

Phase-3 SmartPixels L1 tracking: the `L1SmartPixelsTrackProducer` (SmartPixels-modified
TTTrack collections: passthrough / regression-smeared / future fitted variants), the
customisation entry points for producing and injecting those collections, and the
standalone STEP1/STEP2 drivers. The companion NanoAOD tables live in
`DPGAnalysis/Phase3SmartPixelsNanoAOD`.

## Setup (CMSSW_20_1_X / master cycle)

correctionlib is a proper CMSSW external now (no scram-venv/pip or
SCRAM_CXX_STANDARD/USER_CXXFLAGS hacks); `plugins/BuildFile.xml` declares
`<use name="correctionlib"/>` directly.

```bash
# SCRAM_ARCH=el9_amd64_gcc14
scram project CMSSW CMSSW_20_1_0_pre1
cd CMSSW_20_1_0_pre1/src
cmsenv
git cms-init --upstream-only -q -y
git remote add njmanganelli-fork https://github.com/NJManganelli/cmssw.git
git fetch njmanganelli-fork smartpixels-phase3
# smartpixels-phase3 = master + cms-sw/cmssw#51503 (L1TK community update)
#   + l1nano-central (general L1Nano extensions)
#   + L1Trigger/Phase3SmartPixels + DPGAnalysis/Phase3SmartPixelsNanoAOD
git checkout smartpixels-phase3

# resolution corrections used by the correctionlib* emulator modes
xrdcp root://cmseos.fnal.gov//eos/uscms/store/group/lpcsmartpixelscms/toysim/resolution/spixel_smear_all_configs_barrel_CalV1_v2p1_compound.json L1Trigger/Phase3SmartPixels/data/

scram b -j8
```

## Workflows

Two streamlined workflows are provided by
`L1Trigger/Phase3SmartPixels/python/customizeSmartPixels_cff.py`
(see `DPGAnalysis/Phase3SmartPixelsNanoAOD` for the nano flavors and full cmsDriver recipes):

- **Coexist (WF1)**: one L1Nano file containing BOTH the standard L1Tracks tables and
  SmartPixels L1Tracks tables for one or more variants — in-file track-to-track comparisons.
- **Co-opt (WF2, "pure view")**: `injectSmartPixelsTrackProducer` rewires every downstream
  consumer of the standard tracklet tracks (GTT input, vertex finder, correlator layers,
  track jets/HT/MET, GT, muon matching, ...) to ONE SmartPixels variant; any nano flavor then
  reflects that single interpretation. Comparisons are file-to-file against a standard-tracks
  baseline job.

## Truth postures (`truthSource`)

The TT truth-association maps are Ref/Ptr-keyed to specific track/cluster/stub
collections; a stale map applied to remade tracks silently passes tracks through
the regression/refit producers unmodified (study-breaking). Both
`smartPixelsCoexist` and `smartPixelsCoopt` therefore take a closed
`truthSource` vocabulary of three postures (see `doc/PostureGapStudy.md`):

- **`inJob`** (default) — run the TT associators in-job (unscheduled). Valid only
  when the job also runs DIGI (`DIGI:pdigi_valid,L1TrackTrigger,...`) or the input
  retained `mix:Tracker` `PixelDigiSimLink`s. Fresh new-layout tracks with real
  covariance, but PU is destroyed by re-digitizing signal-only `g4SimHits`.
  Use for: no-PU RelVals, or any sample where re-digitization is acceptable.
- **`fromFile`** — remove all in-process associators; read STEP1-style maps
  straight from the file. Real PU everywhere, but the file's tracks are the old
  (pre-PR#51503) layout so `helixCovMat` is all-zero. Use for: PU studies against
  the file's *own* downstream objects (file L1 jets, etc.), and old-release track
  studies. `digiRefit` here **must** use `seedCovMode="parametrized"`.
- **`fromFileStubs`** (posture C) — rebuild NEW-layout tracks from the file's
  persisted stub tier (`ProducerDTC` → tracklet emulator(s) → re-run track
  associator vs the file's cluster/stub maps); remove only the cluster/stub
  associators, never DIGI. Real PU **and** real per-track covariance, at
  ~2-3 s/event with no pileup re-mixing. `extendedTracks=True` (default) also
  rebuilds the displaced chain. Use for: the primary PU development posture and
  any trackCov-seeded refit production.

`seedCovMode="trackCov"` validity (digiRefit default):

| truthSource | tracks | `trackCov` valid? |
|---|---|---|
| `inJob` | fresh, new layout | **yes** (no PU) |
| `fromFile` | file, old layout | **no** — zero cov; runtime guard `SmartPixelsSeedCovMissing` fires. Use `parametrized`. |
| `fromFileStubs` | fresh, new layout | **yes** (with PU) |

The `fromFile` + `trackCov` combination is deliberately *not* special-cased at
config time; it fails loudly at runtime via the `SmartPixelsSeedCovMissing`
`endStream` guard, by design. Use `fromFileStubs` for trackCov on PU files.

## Refit modes (tier model)

Beyond the passthrough/regression variants, `customizeSmartPixels_cff.py` reserves a
two-tier refit model that develops the OT-L1Track ⊕ IT-pixel refit before the true
SmartPixels chain exists:

- **`digiRefit` (Tier 2, interim)** — project OT L1Tracks into the inner barrel (BPIX)
  pixel layers and refit against the **real** `simSiPixelDigis` digis (classified via
  `PixelDigiSimLink`), **synthesizing only the angle information** (PixelAV response). It is
  an *active-layer* mode (reuses the activeSP `AAII` encoding to pick the refit layer set,
  e.g. `"1100"`) and is *truth-required* (needs pixel digis + `PixelDigiSimLink` +
  TrackingParticles; posture-B/`inJob` or file-present IT products). The producer is
  implemented (Phase 2); the config surface is `DIGIREFIT_DEFAULTS` + the `digiRefitConfig=`
  kwarg on `smartPixelsCoexist`/`smartPixelsCoopt`, and a `digiRefit` variant additionally
  emits the refit sidecar product described below.
- **`refit` (Tier 3, reserved)** — the true system: the producer will ingest OT L1Tracks
  plus a real `SmartTracklet` collection from an `L1SmartTracksFinder`. The mode name is
  accepted by the vocabulary but building a variant with it raises `NotImplementedError`
  (reserved); use `digiRefit` for the interim.

FPGA-fidelity handles live in `DIGIREFIT_DEFAULTS` from day one (float impl, every
truncation switchable to chart resolution-vs-fidelity curves): `useAngles`
(`none`/`alpha`/`alphaBeta`), `maxHitsPerWindow` (combinatorics truncation), `maxKFUpdates`
(Kalman-update cap), and `gainMode` (`full`/`lut` table-driven placeholder).

### digiRefit sidecar product and RNG scheme

`digiRefit` emits a `smartpixels::SmartPixelsRefitSidecar` next to the refit track
collection (same module, same instance label; prompt and extended each emit their own).
It is the single persistent home for every NEW SmartPixels fact absent from the OT-only
`TTTrack`: one `SmartPixelsRefitHitInfo` per attempted layer crossing (residuals,
synthesized angles + sigmas, KF pulls, per-crossing chi2 increments, window multiplicity,
truth-only simlink class + parent angles) plus one `SmartPixelsRefitTrackInfo` per track
(status bits, hit/update counts, chi2 totals). The authoritative contract is
`doc/RefitSidecarSpec.md` (spec v0); the quantizers and the 16-bit "compact" transmitted
word live header-only in `interface/SmartPixelsTransmittedSubset.h` (the single source of
truth ngtagger-train mirrors). A **1:1 output-sync invariant** (one output track and one
sidecar row per input track, in order) is asserted before `put()` and throws
`cms::Exception("SmartPixelsSyncBroken")` on any size mismatch — so PF/Puppi track refs,
DispVertex indices, and the truth table stay valid row-wise against the refit collection.

Randomness (angle synthesis) uses a **local** `CLHEP::MixMaxRng` constructed per event and
seeded from `hash(module label, run, lumi, event)` (FNV-1a). This replaces the old
per-label `RandomNumberGeneratorService` stream engine: outputs are now
event-order-independent and split-job invariant (a physics event yields the same refit
result no matter which job/stream/position processes it), which is what bitwise-reproducible
training productions require. No `RandomNumberGeneratorService` is wired for `digiRefit`.

Migration of the jet-tagging workflow (STEP1/STEP2 with the NGJet tagger) to
20_1/master is in progress; until that lands, the 15_1_0_pre1 instructions below
remain the reference for jet-tagging studies.

## 15_1_0_pre1 instructions with jet tagging for v2p1 of SmartPixels studies
Follow nominal instructions from installer script
```
#!/bin/bash

CMSSW_VERSION=CMSSW_15_1_0_pre1
CMSSW_FORK=CMS-L1T-Jet-Tagging:for-DPS

scram p CMSSW ${CMSSW_VERSION}
cd ${CMSSW_VERSION}/src
[ "$?" != "0" ] && { echo "Unable to set up CMSSW" >&2 ; exit 1; }
eval $(scram runtime -sh)
git cms-init  --upstream-only -q -y
echo "git cms-checkout-topic -u ${CMSSW_FORK}"
git cms-checkout-topic -u ${CMSSW_FORK}
echo "git remote add fork https://github.com/${CMSSW_FORK%%:*}/cmssw.git -t ${CMSSW_FORK##*:} -f"
git remote add fork https://github.com/${CMSSW_FORK%%:*}/cmssw.git -t ${CMSSW_FORK##*:} -f 2>&1 | grep -v 'new tag.*CMSSW'

git clone --quiet https://github.com/cms-hls4ml/hls4mlEmulatorExtras.git && \
  cd hls4mlEmulatorExtras &&
  git checkout -b v1.1.3 tags/v1.1.3
make
make install
cd ..
git clone --quiet https://github.com/Xilinx/HLS_arbitrary_Precision_Types.git hls

git clone https://github.com/CMS-L1T-Jet-Tagging/FastPUPPI.git -b 15_1_0/L1TSC4NGJetTagger
cd FastPUPPI
git reset --hard 8739780aeb4d7d761fbdb07eed1a77798cd62cfd
cd ../../..
```

Then augment with following (rough) instructions, take care with order, it's all very fragile
```
In src directory:
xrdcp root://cmseos.fnal.gov//store/group/lpcsmartpixelscms/jettagging/L1TSC4NGJetModels/L1TSC4NGJetModels.tar.gz .
xrdcp root://cmseos.fnal.gov//eos/uscms/store/group/lpcsmartpixelscms/toysim/resolution/spixel_smear_all_configs_barrel_CalV1_v2p1_compound.json L1Trigger/TrackFindingTracklet/data/
scram-venv cmsenv
cmsenv
python3 -m pip install --no-binary=correctionlib 'correctionlib>2.6' xxhash
git clone https://github.com/cms-hls4ml/L1TSC4NGJetModel.git
git clone https://github.com/CMS-L1T-Jet-Tagging/hls4ml-jettagger.git
mv hls4ml-jettagger L1TSC4NGJetModel
cd L1TSC4NGJetModel
make install
cd ..
# our customization to try and add in the SmartPixels for input production (step 1) too...
git cms-addpkg L1Trigger/TrackFindingTracklet
scram b clean && scram b -j8
# Make sure the runJetNtuple.py works, if needed point to the built modely by adding
# Patch the tagger path since this one isn't in cms-externals for this version of CMSSW, but we built it manually following https://codimd.web.cern.ch/pB3K4fFiSrmblUHFAMYoxA#
# process.l1tSC4NGJetProducer.l1tSC4NGJetModelPath = cms.string(os.environ['CMSSW_BASE']+"/src/L1TSC4NGJetModel/L1TSC4NGJetModel_v0")
git remote add njmanganelli-fork git@github.com:NJManganelli/cmssw.git
git fetch njmanganelli-fork CMS-L1T-Jet-Tagging_for-DPS_cmssw_15_1_0_pre1_v2p1_smartpixelstrackproducer
git checkout CMS-L1T-Jet-Tagging_for-DPS_cmssw_15_1_0_pre1_v2p1_smartpixelstrackproducer
export SCRAM_CXX_STANDARD=20; export USER_CXXFLAGS="$(correction config --cflags)"; export USER_CXXFLAGS=${USER_CXXFLAGS:s/17/20/};export USER_LDFLAGS="$(correction con\fig --ldflags --rpath)"; env | grep 'USER_CXXFLAGS\|USER_LDFLAGS\|SCRAM_CXX_STANDARD';
scram b
# this base-script is for the regular runJetNtuple.py
cmsRun FastPUPPI/NtupleProducer/python/runJetNTuple.py
STEP1 example:
# smartpixels augmented inputs production with ALL variations currently supported running concurrently so that the faster runJetNtuple can reference them (note the module labels having "W" in the name to separate the regular/extended producer from the version, and the replacement of "0" with "I" for inactive and "1" with "A" for active smartpix layers in correctionlib* modes
cmsRun -n 4 L1Trigger/TrackFindingTracklet/python/runInputs151X_withAllSmartPixelsCollections.py inputFiles=/store/mc/Phase2Spring24DIGIRECOMiniAOD/TT_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_AllTP_140X_mcRun4_realistic_v4-v1/2560000/e5b69706-1a56-4df6-b18a-f79a7e4b0001.root outputFile=test.root
STEP2 example:
cmsRun -n 4 L1Trigger/TrackFindingTracklet/python/runJetNTuple_withSmartPixelsInjected.py emulatorMode=correctionlibRegression activeSP=1100
```

## OLD 14_2_0_pre1 instructions
### Ensure good correctionlib version
You must have correctionlib version > 2.2, this has been tested (as of writing) with 2.7 installed in-place of the default available with CMSSW 14. Ensure that an older version is not being picked up or linked against.

### Setup CMSSW
```bash
# Get CMSSW
echo $SCRAM_ARCH
cmsrel CMSSW_14_2_0_pre1
cd CMSSW_14_2_0_pre1/src
cmsenv
# Activate venv
scram-venv cmsenv
cmsenv

# install correctionlib
python3 -m pip install --no-binary=correctionlib 'correctionlib>2.6'

# Setup packages
git cms-init
git cms-addpkg DQMOffline/L1Trigger DataFormats/L1TCalorimeterPhase2 DataFormats/L1TCorrelator DataFormats/L1TParticleFlow DataFormats/L1TrackTrigger DataFormats/L1Trigger DataFormats/StdDictionaries L1Trigger/Configuration L1Trigger/DemonstratorTools L1Trigger/L1TNtuples L1Trigger/L1TTrackMatch L1Trigger/Phase2L1ParticleFlow L1Trigger/Phase2L1Taus L1Trigger/TrackFindingTracklet L1Trigger/TrackTrigger L1Trigger/VertexFinder PhysicsTools/NanoAOD SimTracker/TrackTriggerAssociation
git remote add njmanganelli-fork git@github.com:NJManganelli/cmssw.git
git fetch njmanganelli-fork cmssw_14_2_0_pre1_smartpixelstrackproducer
git checkout cmssw_14_2_0_pre1_smartpixelstrackproducer

# Add L1Nano repo (since merged into CMSSW, but code and instructions need to be adapted to central variation)
cd $CMSSW_BASE/src
git clone git@github.com:cms-l1-dpg/Phase2-L1Nano.git PhysicsTools/L1Nano
cd PhysicsTools/L1Nano
git remote add njmanganelli-fork git@github.com:NJManganelli/Phase2-L1Nano
git fetch njmanganelli-fork vX_1420pre1_gtttracks
get checkout vX_1420pre1_gtttracks
cd -

# compile with C++20, override correctionlib in case it returns -std=c++17, this is a hack until the necessary flags are in a BuildFile.xml
cd $CMSSW_BASE/src
export SCRAM_CXX_STANDARD=20; export USER_CXXFLAGS="$(correction config --cflags)"; export USER_CXXFLAGS=${USER_CXXFLAGS:s/17/20/};export USER_LDFLAGS="$(correction con\fig --ldflags --rpath)"; env | grep 'USER_CXXFLAGS\|USER_LDFLAGS\|SCRAM_CXX_STANDARD'; scram b clean && scram b -j8
```