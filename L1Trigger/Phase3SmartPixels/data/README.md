# Data files (not tracked in git)

The correctionlib* emulator modes of `L1SmartPixelsTrackProducer` need the SmartPixels
resolution-correction set copied into this directory:

```bash
xrdcp root://cmseos.fnal.gov//eos/uscms/store/group/lpcsmartpixelscms/toysim/resolution/spixel_smear_all_configs_barrel_CalV1_v2p1_compound.json L1Trigger/Phase3SmartPixels/data/
```

## cluster_cnn/ (tracked)

The `angleSource = "cnn"` mode of `SmartPixelsRecHitProducer` (opt-in via
`digiRefitConfig={'angleSource': 'cnn'}`) reads these three files:
- `spix_cnn_barrel_v1.onnx`: the TBPX model. It is the `base` variant of
  ngtagger-train `eval_spixel_angles/cluster_cnn`, exported by `export_onnx.py`.
- `spix_cnn_disc_v1.onnx`: the TFPX/TEPX model (the `disc` variant).
- `spix_cnn_angle_sigma_v1.json`: the static sigAlpha/sigBeta tables, written by `sigma_table.py`.

The scalar-input lists in `customizeSmartPixels_cff.SPIX_CNN_ANGLE_SOURCE` must match
each model's `norm.json` `cols`.
