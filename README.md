# Quantized Inertia: Framework-Independent Audit and Conditional TCGS Bridge

Henry Arellano-Peña · Nuevo Estándar Biotropical NEBIOT S.A.S., Bogotá, Colombia
<harellano@unal.edu.co>

Two companion manuscripts, released together under CC BY 4.0.

---

## Contents

| Path | What it is |
|---|---|
| `manuscripts/main.tex` | **Quantized Inertia under Covariant and Quantum-Field-Theoretic Audit: Proper Acceleration, Rindler Stress, Particle-Mass Claims, and Falsifiable Tests** — framework-independent, uses no TCGS–SEQUENTION content |
| `manuscripts/tcgs_bridge_companion.tex` | **Non-Exhaustive Common-Source Inertial Bridges in TCGS–SEQUENTION: Identifiable Low-Rank Images, Anchored Coordinates, and Lossy Cross-Channel Tests** — conditional companion |
| `manuscripts/*.bib` | BibTeX databases (kept beside the `.tex` so `\bibliography{}` resolves) |
| `figures/` | Generated figures |
| `scripts/qi_scalar_response.py` | Regenerates Figure 1 of the audit |
| `CHANGELOG.md` | Version history |

The two papers are deliberately firewalled. The audit establishes its
quantum-field-theoretic, covariance, algebraic, laboratory, and astrophysical
conclusions **without** using TCGS–SEQUENTION, and the companion does not use
their correctness as evidence for Counterspace. Neither is permitted to borrow
evidentiary status from the other.

---

## What the audit concludes

Quantized inertia is a **concrete and falsifiable but presently incomplete**
modified-inertia programme. Firm results, each targeting a specific derivation
rather than the programme as a whole:

- The Minkowski vacuum restricted to a Rindler wedge is KMS with respect to
  boosts, but its local renormalized stress tensor remains zero — a Rindler
  decomposition alone supplies no pressure imbalance.
- The Fulling–Rindler vacuum has negative comoving energy density
  `−1/(480π²ξ⁴)` and **vanishing** momentum flux; it supplies no directed thrust,
  and is not the state produced merely by accelerating through the Minkowski
  vacuum.
- The Unruh relation uses proper acceleration; QI's galactic applications use
  acceleration relative to surrounding matter. Without a covariant bridge the
  laboratory and galactic branches are not two tests of one law.
- The scalar law is **implicit** in `a`. Closing it on the realized acceleration
  gives `a = F/m + a₀` and `mᵢ/m = F/(F + ma₀) ∈ [0,1)`, so positivity and the
  minimum acceleration are theorems, not appended admissibility rules. The
  frequently repeated "negative inertial mass" objection applies only to the
  explicit kinematic reading and is not pressed here.
- The primary galactic derivation applies that closure to a **decomposed**
  acceleration `|A| = |a| + |a′|`. Both published relations — the acceleration
  floor and the Tully–Fisher relation — are specializations of one master
  equation. The physical work is done by an undeclared decomposition rule, whose
  summation domain is never specified; that is also what shields the floor from
  every laboratory bound.
- The assignment `a′ = a₀` at a galaxy's edge is asserted, not derived, and
  carries a sharp untested consequence: internal small-scale accelerations must
  **rise** in the outer disc. This is checkable against resolved kinematics and
  is **not shared with MOND** — the strongest galactic test the programme offers.
- Hamiltonian variance is not automatically extractable energy; `ΔE Δx ~ ℏc` has
  limited spectral jurisdiction; the gravity derivation imports `G` through the
  Planck mass; `ℏ/(cλ_C) = m/2π`, while the reduced wavelength returns the mass
  by definition.
- The capacitor coefficient `α = π²cmₑ/(120e) = 1.4019×10⁻⁴ N A⁻¹` is
  reconstructible only under the source's half-final-speed convention; the
  final-speed convention gives a value four times smaller.
- The specific QI wide-binary peak prediction is rejected by the Gaia DR3
  distributional test. The broader wide-binary gravity dispute remains open in
  **both** directions.

## What the companion establishes

A source–Shadow framework acquires no empirical content by factoring carrier
responses through an unrestricted latent operator — any finite family of
unrelated kernels admits an automatic block-diagonal realization. The companion
constructs a bridge that is genuinely restrictive: capacity counted on the
**identifiable image** rather than the nominal rank variety, a preregistered
compression ceiling and parameter-to-data ceiling, one global coordinate gauge
with per-channel gauges forbidden, metric-anchored pseudoinverses, and a powered
held-out test. A reproducible three-channel synthetic example accepts a true
rank-two bridge at the declared 5% size and rejects independent kernels with unit
power; a near-saturated construction retaining 98% of unrestricted capacity is
rejected by the compression ceiling.

Success supports **that bridge only**. It does not establish Counterspace, the
Extrinsic Constitutive Law, or non-ontic time.

---

## Building

Both documents are standard `article` class and build with pdfLaTeX + BibTeX:

```bash
cd manuscripts
pdflatex main && bibtex main && pdflatex main && pdflatex main
pdflatex tcgs_bridge_companion && bibtex tcgs_bridge_companion \
  && pdflatex tcgs_bridge_companion && pdflatex tcgs_bridge_companion
```

`\graphicspath{{figures/}}` is set in both preambles, so either build from the
repository root or copy `figures/` beside the sources. On Overleaf, upload the
tree as-is.

Required packages: `geometry, fontenc, inputenc, lmodern, microtype, setspace,
amsmath, amssymb, amsthm, mathtools, bm, mathrsfs, siunitx, graphicx, booktabs,
array, tabularx, longtable, float, caption, enumitem, tcolorbox, natbib,
hyperref, cleveref` (companion also uses `tikz`).

### Figures

`figures/qi_scalar_response.pdf` is generated by
`python3 scripts/qi_scalar_response.py` (needs `numpy`, `matplotlib`).

**Not yet in this repository** — supply from your local working copy:

- `figures/qi_capacitor_scaling.pdf`
- `figures/bridge_synthetic_residuals.pdf`

Also not yet included: `code/bridge_synthetic_demo.py` and its JSON result file,
which the companion's Appendix C and Data Availability statement both promise.
Adding them is a prerequisite for the reproducibility claim.

---

## Before submission

- [ ] Mint persistent identifiers (Zenodo or Preprints.org) for both manuscripts
      so they cite each other resolvably. The companion's `\cite{ArellanoQIAudit2026}`
      is load-bearing for its branch-separation argument. See the flagged block at
      the top of the companion `.bib`.
- [ ] Add the two missing figures and the synthetic-demo script above.
- [ ] Bump `\date{}` and the version note in the companion `.bib`.
- [ ] Verify §2.3 of the audit against McCulloch (2012) eqs. 2–9 one final time —
      it is a close reconstruction of a primary derivation and carries the paper's
      sharpest claim.

## License

CC BY 4.0 — see [`LICENSE`](LICENSE). Externally cited works retain their own
licenses.
