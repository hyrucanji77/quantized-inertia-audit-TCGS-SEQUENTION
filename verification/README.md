# Numerical verification and provenance

This directory contains numerical evidence for the V8 computational release.
The primary executable programs are in [`code/`](../code/) and
[`scripts/`](../scripts/). No manuscript source or source-diff verification is
needed to run them.

## Included records

| File or directory | Purpose |
|---|---|
| [`three_channel_statistical_preservation.json`](three_channel_statistical_preservation.json) | Confirms preservation of the original three-channel configuration, truth, thresholds, and per-replicate statistics while correcting capacity/data metadata |
| [`bridge_four_channel_verification.json`](bridge_four_channel_verification.json) | Independent vectorized least-squares checks and integrity checks on recorded four-channel results |
| [`bridge_four_channel_generating_v8.py`](bridge_four_channel_generating_v8.py) | Exact original source identified by the primary four-channel design and result records |
| [`bridge_four_channel_loader_amendment.json`](bridge_four_channel_loader_amendment.json) | Source hashes, the implementation-only loader correction, tolerance guard checks, and comparison with the complete rerun |
| [`four_channel_loader_rerun/results/bridge_four_channel_design.json`](four_channel_loader_rerun/results/bridge_four_channel_design.json) | Frozen design used by the separate rerun |
| [`four_channel_loader_rerun/results/bridge_four_channel_results.json`](four_channel_loader_rerun/results/bridge_four_channel_results.json) | Complete rerun under the active executable, including current and original frozen-script hashes |
| [`four_channel_loader_rerun/figures/bridge_four_channel_results.pdf`](four_channel_loader_rerun/figures/bridge_four_channel_results.pdf) | Regenerated figure from that rerun |

## Original generating code and active executable

The primary four-channel result file preserves the hash of the exact program
that generated it. That source is archived here without changing its bytes:

```text
verification/bridge_four_channel_generating_v8.py
SHA-256: 92c9bdc03b763f630ab2859036a71aca1bd890d30eb754230efe6d3636313a1c
```

Use the active executable for reproduction:

```text
code/bridge_four_channel_demo.py
SHA-256: f5b886df757aec5f67d1aa34444c4ff94c1b94e6ecb0977c8b81c25df441b76b
```

The active version checks the discrete configuration exactly, validates regenerated
fixed arrays with relative and absolute tolerances of `2e-12`, and then uses the
saved arrays. Structural and least-squares diagnostics are recomputed under their
declared tolerance assertions instead of being compared for bit identity.
This avoids rejecting a mathematically equivalent design because of harmless
QR/SVD roundoff.

The original frozen design, primary result file, and primary figure remain
unchanged. The separate 7,500-replicate rerun reproduced every recorded statistic
and every alternative readable block exactly in the recorded environment. Its
output distinguishes the current executable hash from the hash carried by the
original frozen specification. Metadata such as run timestamps is not expected
to be identical on a later rerun.

## Checks and interpretation

The four-channel portfolio operator has rank nine, and its image on the regular
rank-two tangent family has dimension eight. The recorded alternatives preserve
rank two, the fixed input/output ranges, and squared Frobenius norm `1.25` per
readable block. Independent vectorized least-squares calculations check the
implemented pooled and per-channel estimators.

For the retained three-channel experiment, all per-replicate statistic arrays
and reported thresholds and rejection summaries are unchanged. One diagnostic
Gaussian SVD differs by at most approximately `8.53e-14` in singular values,
without changing its rank. Cross-platform verification should therefore use
appropriate numerical tolerances rather than require identical floating-point
diagnostic bytes.

These records verify the specified calculations and implementation. They do not
establish a physical inertial mechanism, uniform size over an entire composite
null, local power for the four-channel design, or completion of every physical
bridge-protocol condition. The [main README](../README.md) states the model and
capacity distinctions governing the numerical results.

The author-generated code and verification records are released under CC BY 4.0;
see [`LICENSE`](../LICENSE).
