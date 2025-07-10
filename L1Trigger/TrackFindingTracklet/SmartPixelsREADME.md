
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