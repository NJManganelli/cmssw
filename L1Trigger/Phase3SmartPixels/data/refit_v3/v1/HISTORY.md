# Payload set v1: byte history

The set id stays v1 across the change below: no value changed and nothing produced had consumed the set.
Results that depend on the set name the skip_price.json checksum, so the two byte versions stay distinguishable.

| export (UTC) | commit | skip_price.json sha256 | md5 | content |
|---|---|---|---|---|
| 2026-10-20 12:05:56 | ed9b8522eab | 99a036cd2ce74febebcfab4b703a0fc3235073d3c3e428213b29ef33e4d139c1 | 258ad4347e185090d491b3589864a170 | skip_price |
| 2026-10-20 20:41:41 | 7ea72fd5b63 | 07eec3fc87310b50ec6a82d1f8075d85531aa2a76eafe6c8541a7b17e23c7599 | d1c70615960d77cd18ac692939493468 | skip_price, skip_price_level |

- Why the bytes changed: the file gained the diagnostic correction `skip_price_level` (the level of the cut
  file's fallback chain that gives each skip price), read only for the closure's per-path counts. The `skip_price`
  correction is the identical JSON object in both; every other file of the set is byte-identical.
- An uncommitted export at 20:29:44 (used by the side-lane job v3-pregate-1) came from the same table code. The
  export is byte-deterministic (a re-export reproduces 07eec3fc...), so it carried the 20:41 bytes.
- The v3 dev jobs before 20:29 used the first bytes. The v3 job of file 4 run with each version is bitwise equal on
  every common nano branch.
