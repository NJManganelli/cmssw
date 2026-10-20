# Seedless refit v3: re-calibration procedure

How a new payload set `data/refit_v3/v<N>/` is made, for example when the real SmartPixels sensor's
angles replace the emulator's placeholder angles. Each step is one command; the reference
implementation is ngtagger-train `eval_refitq/replay/`.

1. **The new angle source** goes into `SmartPixelsRecHitProducer` (the only place angles are made) and
   a production (canon campaign) is made with it.
2. **The angle tables:** refit the barrel (global key) and disc (module-frame key) angle-error tables
   with the locked P0.3 recipe on the new production (`eval_refitq/p03_angles`).
3. **The skip price:** refit the cut tables with the replay's cut fit, all four dims sets
   (`gate|xy`, `xya`, `xyb`, `xyab`) (`eval_refitq/replay/run_c1b.py` J4). The seed model
   (`replay/seedmodel.py`) and the material scale (`run_c1b.py` J3) are position-only fits: refit them
   only when the positions, the geometry or the OT track finding change.
4. **Physics checks, required:** re-run the disc calibration's angle pulls per key bin and its cut
   coverage on the replay with the new tables, and register them. Their registry ids go into the
   manifest (`--physics-checks`); no set ships without them.
5. **Export:** the replay agent exports the specification (`replay/export_v3.py`), then
   ```
   python eval_refitq/v3cmssw/export_cmssw.py --spec <spec> --set-id v<N> --physics-checks <ids> \
       --out L1Trigger/Phase3SmartPixels/data/refit_v3/v<N> \
       --unit-ref L1Trigger/Phase3SmartPixels/test/refit_v3_unit_reference.json
   ```
   The converter tabulates every lookup through the replay's own code at every index along each axis
   and checks the written correctionlib files against itself before writing the manifest.
6. **Closure:** `testSmartPixelsRefitV3` (unit) and `refitV3NanoClosure` (the engine against the
   replay's reference on 10 events per scenario) must pass, then point `seedlessRefitV3.payloadManifest`
   at the new set.
