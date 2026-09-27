# Data files (not tracked in git)

The correctionlib* emulator modes of `L1SmartPixelsTrackProducer` need the SmartPixels
resolution-correction set copied into this directory:

```bash
xrdcp root://cmseos.fnal.gov//eos/uscms/store/group/lpcsmartpixelscms/toysim/resolution/spixel_smear_all_configs_barrel_CalV1_v2p1_compound.json L1Trigger/Phase3SmartPixels/data/
```

## smartpixels_emu/

The `angleSource = "emu"` mode of `SmartPixelsRecHitProducer` (the SmartPixels-output
emulator; opt-in via `digiRefitConfig={'angleSource': 'emu'}`) reads six files, all tracked
and all made by ngtagger-train `eval_spixel_angles/cluster_cnn` (`emu_*.py`):
- `spix_emu_cnn_barrel_v2.onnx`, `spix_emu_cnn_disc_v2.onnx`: the tiny cluster CNNs (fast-study
  models `barrel_tiny`, `disc_tiny`), exported by `emu_export_onnx.py cnn`. The scalar-input
  lists in `customizeSmartPixels_cff.SPIX_EMU_ANGLE_SOURCE` must match each model's
  `norm.json` `cols`.
- `spix_emu_a1_barrel_v2.onnx`, `spix_emu_a1_disc_v2.onnx`: the A1 class-C residual models
  (trained by `emu_train_a1.py` on class T of the 125 training events), with the draw inside
  the graph, exported by `emu_export_onnx.py a1`.
- `spix_emu_library_v2.root` (1.5 MB): the rescaled SmartPixels Conv1D_Full-2bit dump (tree
  `dump`), written by `emu_build_library.py`.
- `spix_emu_sigma_v3.json`: the static sigAlpha/sigBeta tables, by surface code, cluster size
  along the angle's axis (sizeX for alpha, sizeY for beta; the size is taken as sensor-reported)
  and |reported cot|, from the emulated class-T and class-C routes of the training events
  weighted by the class fractions of each (surface, size) bin, written by `emu_sigma_table.py
  build ... own`. The producer takes the inputs by name (surface, sizeX, sizeY, absCot).
