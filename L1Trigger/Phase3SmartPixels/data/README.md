# Data files (not tracked in git)

The correctionlib* emulator modes of `L1SmartPixelsTrackProducer` need the SmartPixels
resolution-correction set copied into this directory:

```bash
xrdcp root://cmseos.fnal.gov//eos/uscms/store/group/lpcsmartpixelscms/toysim/resolution/spixel_smear_all_configs_barrel_CalV1_v2p1_compound.json L1Trigger/Phase3SmartPixels/data/
```

## spix_angle_response_Conv1D_Full-2bit_v6.json, _v5.json

The PixelAV angle-response payload of `angleSource = "payload"` (`seedlessRefitConfig['pixelavAngleSet']`;
v6 is the default of `test/makeSpixConfig.py`): bias and sigma of the SmartPixels Conv1D_Full-2bit
network in bins of the true (cotAlpha, cotBeta), CMS frame, written by ngtagger-train
`eval_spixel_angles/extract_pixelav_angle_payload.py --version vN --binning vN`. v5 is the first
version in PHYSICAL units: the eval dump stores normalised labels (physical = 1.8827 x stored for CMS
alpha, 6.5289 x stored for CMS beta), which payloads up to v4 used as they were. v5 has 6 cotBeta
bins (|beta| 0-63 deg is one bin). v6 has 24 cotBeta bins, symmetric in the angle (edges at 0, 15,
25, 35, 45, 50, ..., 80 deg, then the v5 outer edge |cotBeta| = 7.18), grown from the dump so that
every (cotAlpha, cotBeta) cell holds >= 100 dump rows; its top cotAlpha bin is [0.75, 1.51] (v5
split it at 1.13, but the dump ends at 1.40).

## smartpixels_emu/

The `angleSource = "emu"` mode of `SmartPixelsRecHitProducer` (the SmartPixels-output
emulator; opt-in via `seedlessRefitConfig={'angleSource': 'emu'}`) reads six files, all tracked
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
- `spix_emu_sigma_v2.json`: the static sigAlpha/sigBeta tables, by surface code and |reported
  cot|, from the emulated class-T and class-C routes of the training events, written by
  `emu_sigma_table.py build ... none`. The producer takes the inputs by name (surface, sizeX,
  sizeY, absCot); a size-keyed table (`own`/`both`, v3) was studied and not adopted.
