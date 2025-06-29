
### Ensure good correctionlib version
You must have correctionlib version > 2.2, this has been tested (as of writing) with 2.5 as available in CMSSW_14_2_0_pre1. Ensure that an older version is not being picked up or linked against.

### Setup CMSSW

### Setup venv (1-time setup)
```bash
scram-venv cmsenv
cmsenv
```


```bash
#compile with C++20, override correctionlib in case it returns -std=c++17
export SCRAM_CXX_STANDARD=20 export USER_CXXFLAGS="$(correction config --cflags)" export USER_CXXFLAGS=${USER_CXXFLAGS:s/17/20/} export USER_LDFLAGS="$(correction config --ldflags --rpath)"; print $SCRAM_CXX_STANDARD "\n" ${USER_CXXFLAGS} "\n" $USER_LDFLAGS; scram b clean && scram b -j8
```