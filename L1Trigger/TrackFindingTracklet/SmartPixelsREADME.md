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