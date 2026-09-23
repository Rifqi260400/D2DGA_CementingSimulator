# Archived: the evidence behind the pre-2026-09-19 K-GEP-1 numbers

`kgep1_ckpt_synthetic_no_axial_gradient_16x80_w0.2_v1.2.npz` (2026-09-17) is the
run whose `eta_E = 0.9946` appears in `AUDIT_REPORT.md`, `docs/gate_status.md`
and `output/kgep1_results.md`. It predates both `d2dga/runio.py` (recorded
provenance, R-2) and `d2dga/conservation.py` (the BCF25 §IV ledger, REM-16), so:

* its settings are **reconstructed** from the filename plus the current
  `config.py`, not recorded;
* its conservation figure, `2.93e-14`, uses the **older normalisation** (by the
  volume present, not by the annulus volume), so it is not the quantity BCF25
  report and cannot be set beside their ~1e-15;
* its closure table was **loaded from cache**, so the B-1 angle sensitivity was
  never measured for it.

`closure_table_b0a79f5ea2da.npz` is that cached table.

Both are kept here rather than deleted, because they are what the published
figures were computed from. The re-run that replaces them writes to the same
checkpoint name — deliberately, since the documents refer to it — and carries a
recorded config, a metrics file, a conservation ledger and a measured B-1.
