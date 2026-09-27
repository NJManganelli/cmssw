# Data files (not tracked in git)

The correctionlib* emulator modes of `L1SmartPixelsTrackProducer` need the SmartPixels
resolution-correction set copied into this directory:

```bash
xrdcp root://cmseos.fnal.gov//eos/uscms/store/group/lpcsmartpixelscms/toysim/resolution/spixel_smear_all_configs_barrel_CalV1_v2p1_compound.json L1Trigger/Phase3SmartPixels/data/
```

## smartpixels_emu/

The `angleSource = "emu"` mode of `SmartPixelsRecHitProducer` (the SmartPixels-output
emulator; opt-in via `digiRefitConfig={'angleSource': 'emu'}`) reads four files, all made
by ngtagger-train `eval_spixel_angles/cluster_cnn` (`emu_*.py`):
- `spix_emu_cnn_barrel_v1.onnx`, `spix_emu_cnn_disc_v1.onnx` (tracked): the pos-hubsig
  cluster CNNs (models `barrel_pos_hubsig`, `disc_pos_hubsig`), exported by
  `emu_export_onnx.py`. The scalar-input lists in
  `customizeSmartPixels_cff.SPIX_EMU_ANGLE_SOURCE` must match each model's `norm.json` `cols`.
- `spix_emu_sigma_v1.json` (tracked): the static sigAlpha/sigBeta tables, by surface code and
  |reported cot|, written by `emu_sigma_table.py`.
- `spix_emu_library_v1.root` (NOT tracked, too large for git): the rescaled SmartPixels
  Conv1D_Full-2bit dump (tree `dump`) and the class-T donor library of 25 training events
  with its standardisation (trees `donors`, `norm`), written by `emu_build_library.py build`.
  Until it is hosted externally, copy it here from
  `/Volumes/WDMac/smartpixels_scratch/cluster_cnn_v2/emu_cmssw/spix_emu_library_v1.root`.
