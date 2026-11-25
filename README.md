# Quantized Inertia: Framework-Independent Audit and Conditional TCGS Bridge

Henry Arellano-Peña · Nuevo Estándar Biotropical NEBIOT S.A.S., Bogotá, Colombia  
<harellano@unal.edu.co>

Computational materials for two companion studies, released under CC BY 4.0.

## Current computational release: V8 — 25 November 2026

This release supplies the Python programs, frozen computational designs, numerical
results, figures, and numerical verification records accompanying the V8 QI audit
and conditional TCGS–SEQUENTION inertial-bridge companion. It includes the scalar
closure and Cassini-proxy calculations, the original three-channel control, the
partial-readability experiment with its local power curve, and the four-channel
control with capacity conditional on the frozen channel ranges.

**Manuscript version boundary.** The existing files under `manuscripts/`, including
the LaTeX and BibTeX files, belong to the earlier repository release. They remain
unchanged and are **not the manuscript sources for this V8 computational
release**. This update does not add V8 LaTeX sources or Overleaf archives. The
current computational definitions and numerical results are documented below and
in the accompanying machine-readable records; older manuscript and changelog
statements must be read in their historical version context.

The two studies have distinct evidentiary jurisdictions. The QI audit examines
specified physical and mathematical claims independently of TCGS–SEQUENTION.
The companion develops a conditional finite inertial bridge within that
framework. The simulations test declared finite statistical constructions; their
success does not establish a physical inertial mechanism or the governing source
architecture.

## Contents

| Path | Current contents |
|---|---|
| `scripts/qi_scalar_response.py` | Two-panel comparison of the explicit kinematic and implicit realized scalar readings |
| `scripts/qi_capacitor_scaling.py` | Capacitor force-scaling figure for the stated model and coefficient convention |
| `scripts/qi_scalar_cassini_check.py` | Scalar-law identities, residual-dependent excess, conditional Cassini proxies, and numerical rounding checks |
| `code/bridge_synthetic_demo.py` | Original fully readable three-channel control with corrected V8 capacity/data metadata |
| `code/bridge_partial_local_demo.py` | Partial reconstruction, an unread bridge sector, correlated shared-sector tests, and the full local power grid |
| `code/bridge_four_channel_demo.py` | Four-channel control with frozen-range comparison capacity and independent rank-two alternatives inside those ranges |
| [`results/`](results/) | All fixed settings, seeds, maps, thresholds, per-replicate results, and uncertainty summaries |
| [`figures/`](figures/) | Five generated PDF figures; see the [figure index](figures/README.md) |
| [`verification/`](verification/) | Numerical checks, original generating-code provenance, and the separate four-channel rerun |
| [`CHANGELOG.md`](CHANGELOG.md) | Current computational changes followed by preserved historical release entries |

## Bridge comparison and capacity

For frozen channel maps, let `r_in,chi` and `r_out,chi` denote their ranks.
The independent comparison model permits each readable channel kernel to vary
independently inside those same frozen ranges. Its dimension is

```text
D_ind = sum_chi r_in,chi * r_out,chi.
```

This is the denominator of the declared capacity ratio. The raw physical-kernel
count `sum_chi p_chi*q_chi` is diagnostic only. Range restrictions already supplied
by map calibration are not counted as compression supplied by the common bridge.

| Design | Visible common rank-two dimension | Conditional independent dimension | Capacity ratio | Decision at `rho_max = 0.25` |
|---|---:|---:|---:|---|
| Three fully readable channels | 8 | 27 | 0.29630 | Fails the capacity gate |
| Three partially readable channels | 16 | 27 | 0.59259 | Fails the capacity gate |
| Four fully readable channels | 8 | 36 | 0.22222 | Passes the capacity gate |

The declared independent comparison has no additional per-channel rank bound, so
the ratio measures the combined common-bridge and low-rank restrictions. Against
independent rank-two readable blocks, the four-channel denominator is instead
`4*8 = 32`, giving `8/32 = 0.25`. Sharing alone meets the declared ceiling exactly;
the comparison with unrestricted readable blocks supplies the additional margin.

With known rank-three output maps and independently known isotropic noise, the
three-channel control has 270 informative scalar calibration responses out of
540 recorded output coordinates. The four-channel control has 360 out of 720.
Their parameter-to-informative-data ratios are `8/270 = 0.02963` and
`8/360 = 0.02222`. These are observation counts, not design-matrix or
Fisher-information ranks. The four-channel specification fixes `eta_max = 0.1`.

The near-saturated examples use full physical map ranks, so conditional and raw
independent dimensions coincide at 200. The 18-dimensional construction has
visible dimension 196 and kernel dimension 128, giving `196/200 = 0.98`; the
14-dimensional comparison has visible dimension 164 and gives `164/200 = 0.82`.

## Recorded numerical results and their scope

### Original three-channel control

The original Monte Carlo configuration and every per-replicate statistic are
retained. Its calibrated threshold is approximately `413.265`, its independent
null rejection rate is `128/2500 = 0.0512`, and its empirical rejection rate for
the specified independent physical-kernel alternative is `2500/2500`.

That alternative can violate the frozen channel ranges. This is an extreme
residual-exclusion control, and its successful residual behavior does not repair
its failed capacity gate. The code's historical experiment identifier records
the experiment's origin; V8 corrects the capacity/data accounting and gives the
figure a version-neutral title.

### Partial-readability experiment

Three rank-three channel pairs read sectors `{1,2,3}`, `{1,2,4}`, and `{1,2,5}` of
a six-dimensional bridge. The sixth coordinate remains unread. The full
portfolio operator has rank 19 and kernel dimension 17; the generic rank-two
family has nominal dimension 20 and visible dimension 16. A distinct rank-two
bridge with identical channel readouts is explicitly recorded.

The correlated shared-sector statistic has eight independent contrast
coordinates. Its Monte Carlo threshold is `15.3517419027`, and its independent
null rejection rate is `126/2500 = 0.0504`. The complete local grid is retained:

```text
delta = 0, 0.005, 0.010, 0.020, 0.027, 0.030,
        0.045, 0.063, 0.080, 0.100, 0.150, 0.200.
```

Each amplitude uses 2,500 replicates. Alternatives remain inside the frozen
ranges and retain rank two generically. `delta` is a relative RMS response
displacement from the fixed baseline, not distance to the entire compatible
family. At `delta = 0.10`, the empirical rejection rate is `2403/2500 = 0.9612`,
the first tested amplitude meeting the 95% planning target. This is a
grid-specific, ensemble-specific result, not a universal minimum detectable
departure. The design demonstrates partial reconstruction and local
shared-compatibility sensitivity while failing the separate capacity ceiling.

### Four-channel control

Each physical channel has six input and six output coordinates. Fixed rank-three
maps read a common `3 x 3` bridge with singular values `(1, 0.5, 0)`. Each channel
receives 30 calibration and 20 held-out standard-normal probes, with independent
Gaussian output noise of known standard deviation `0.03`. The common estimator
is pooled unconstrained least squares followed by rank-two SVD truncation; it is
not asserted to be the rank-constrained maximum-likelihood estimator.

The alternative uses independently oriented rank-two readable blocks. Every
block remains inside its channel's frozen ranges, has squared Frobenius norm
`1.25`, and is fixed across the replicate's calibration and held-out samples.
The threshold, null evaluation, and alternative evaluation each use 2,500
replicates.

| Quantity | Recorded value |
|---|---:|
| Calibrated 95th-percentile threshold | 537.750510 |
| Null rejections | 118 / 2500 |
| Empirical false-rejection rate | 0.0472 |
| Pointwise Wilson 95% interval | [0.03956, 0.05623] |
| Alternative rejections | 2500 / 2500 |
| Empirical power for the declared ensemble | 1.0000 |
| Pointwise Wilson 95% interval | [0.99847, 1.00000] |
| Median null statistic | 482.814 |
| Median alternative statistic | 86539.754 |

The independent fitted comparison uses unrestricted `3 x 3` blocks inside the
same ranges. Its mean held-out difference `T_ind - T_common` is `+22.0031` under
the common truth and `-86730.1350` under the alternative. These predictive-score
comparisons are descriptive, without an additional calibrated evidence decision.

The size estimate is pointwise at the frozen truth. Power concerns the specified
separated alternative ensemble, and the intervals condition on the realized
Monte Carlo threshold. The local power curve above belongs to the partial
design; a local curve has not been measured for this four-channel design. Passing
its capacity/data ceilings and obtaining these residual results does not establish
every bridge-protocol condition or distinguish a physical bridge from a
calibrated standard-physics common response.

## Scalar QI and Cassini calculations

The scalar checks use the declared positive collinear reduction
`g = a * (1 - a_Q/(a+b))`, with `g = F/m`, `a_Q = 2*c^2/Theta`, and prescribed
residual `b`. The code verifies the exact root, the residual-dependent excess
`Delta = a-g`, the inverse residual relation, and the associated numerical
proxies. At fixed finite `b`, the high-acceleration excess approaches `a_Q` even
though the fractional correction tends to zero. For `b = a_Q`, the scalar map
has the simple-MOND algebraic form at the same input and acceleration scale;
that identity does not establish equivalence of complete theories.

The Cassini comparison rescales a specified persistent sunward perturbation in
the near-circular calculation of [Banik and Kroupa (2020)](https://arxiv.org/abs/2006.06000).
The approximately 2,130 km radial component follows from their Eq. (7),
`delta_r = x_0*(cos(phi)-1)`, at `phi = pi` with `x_0 = 1065 km`. It is distinct
from the approximately 7,000 km total in-plane displacement. The output contains
conditional sensitivity proxies, not a tracking-data analysis, a QI ephemeris
fit, or a universal observational exclusion. Source details and inputs are
recorded in [`results/qi_scalar_cassini_check.json`](results/qi_scalar_cassini_check.json).

## Reproducing the computations

The recorded environment is CPython **3.12.14**, NumPy **2.3.5**, and Matplotlib
**3.10.8**. NumPy and Matplotlib are the only third-party dependencies; the
scalar/Cassini check uses the Python standard library. Install the repository's
[`requirements.txt`](requirements.txt) in a Python 3.12 environment.

Run the following commands **from the repository root**:

```bash
python3 scripts/qi_scalar_response.py
python3 scripts/qi_capacitor_scaling.py
python3 scripts/qi_scalar_cassini_check.py
python3 code/bridge_synthetic_demo.py
python3 code/bridge_partial_local_demo.py --design-only
python3 code/bridge_partial_local_demo.py
python3 code/bridge_four_channel_demo.py --design-only
python3 code/bridge_four_channel_demo.py
```

The supplied `figures/` directory must exist for the scalar-response script.
The bridge scripts also accept `--output-root PATH`; the Cassini script accepts
`--output PATH`. Normal runs overwrite their generated result and figure files.
Use a separate output directory or a separate working copy to retain the
distributed records unchanged.

For example, reproduce the four-channel experiment in a separate directory
using its original frozen arrays:

```bash
mkdir -p reproduction/results
cp results/bridge_four_channel_design.json reproduction/results/
python3 code/bridge_four_channel_demo.py --output-root reproduction --design-only
python3 code/bridge_four_channel_demo.py --output-root reproduction
```

The saved computational specifications precede Monte Carlo generation; they are
not claims of public preregistration. Seeds, draw order, estimators, quantile
rules, and software versions are recorded. Floating-point details may vary across
numerical libraries; the delivered four-channel loader verifies the configuration
and fixed arrays with declared tolerances and uses the saved arrays for the run.

## Numerical provenance

The original four-channel numerical record identifies the exact generating code
preserved as [`verification/bridge_four_channel_generating_v8.py`](verification/bridge_four_channel_generating_v8.py).
The active executable includes a loader correction for harmless numerical
roundoff. Its design and all reported numerical results are unchanged; a separate
complete rerun reproduced every recorded statistic and alternative block.
The [verification index](verification/README.md) explains the source hashes,
checks, and retained rerun evidence.

[`SHA256SUMS.txt`](SHA256SUMS.txt) records the checksums of the files supplied in
this computational release. From the repository root, verify them with
`sha256sum --check SHA256SUMS.txt`. The manifest excludes the unchanged historical
manuscript files and excludes itself.

## License

CC BY 4.0 — see [`LICENSE`](LICENSE). Externally cited works retain their own
licenses.
