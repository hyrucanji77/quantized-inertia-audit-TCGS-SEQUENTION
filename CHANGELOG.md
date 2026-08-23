# Changelog

## v5.0.0 — August 2026

Revision of the V4 pair. The headline change is a **retraction and replacement**:
the frequently made "QI needs an ad hoc rule to avoid negative inertial mass"
objection is withdrawn as a reading artefact, and a sharper, sourced criticism
takes its place.

### `main.tex` — Quantized Inertia audit

**Added — the implicit closure (new §2.2).** The scalar law is implicit in `a`.
Reading E (explicit/kinematic) and Reading R (implicit/realized) are distinguished.
Closing `F = mᵢa` gives `F = m(a − a₀)`, hence `a = F/m + a₀` and
`mᵢ/m = F/(F + ma₀) ∈ [0,1)`. Positivity is automatic; the minimum acceleration is
the `F → 0` limit; `a = 0` is never attained, so the singularity of the scalar
factor is not encountered. New **Principle 1 (Closure declaration)**.

**Added — the acceleration-decomposition postulate (new §2.3).** Reconstructed
from McCulloch (2012), eqs. 2–9. The galactic derivation applies the closure to
`|A| = |a| + |a′|`, giving the master equation

    (|a| + |a′| − a₀) a = GM(|a| + |a′|)/r²

Both published relations are specializations: `a′ = 0` gives the acceleration
floor used for JuMBOs and Proxima; `a′ = a₀` at a galaxy's edge cancels the floor
term and gives `v⁴ = 2GMc²/Θ`. Three consequences: (a) `|A|` is a scalar sum of
magnitudes over an unspecified matter set — not a vector sum, not covariant, not
computable from a system description; (b) the floor is carried by `|A|`, not
`F/m`, which is what keeps it below torsion-balance sensitivity, so the
programme's signature prediction is shielded by an unspecified rule; (c) the
`a′ = a₀` assignment implies rising internal small-scale accelerations in the
outer disc — sharp, untested, and not shared with MOND.

**Added — residual-acceleration test (new §9.5).** Develops (c) above into an
executable programme with the three declarations needed to make it falsifiable
(what `a′` physically is, the aperture of the sum, the predicted radial profile).

**Added — sub-threshold observation (new §9.6, Table 2).** Every astrophysical
target sits at or below `a₀ ≈ 2×10⁻¹⁰ m s⁻²`: Proxima ~7×10⁻¹¹, JuMBO@400 AU
~7×10⁻¹¹, galaxy outskirts 10⁻¹¹–10⁻¹⁰. In that regime predictions are dominated
by the single constant `a₀`, so agreement discriminates weakly. Corollary:
preregister the transition region `|a| ~ a₀`, not more deep-low-acceleration
systems.

**Rewritten — §4.5.** From "Positivity and branch closure" to "Positivity is
closure-dependent, and the residual burden is the acceleration floor." States
plainly that an audit pressing the negative branch is auditing Reading E while
the applications use Reading R. Replaces that objection with three compulsory
declarations: the domain of the acceleration sum, whether the floor is
directional, and the status of the decomposition assumption.

**Updated — §9.2 wide binaries.** Added Chae & Yoon (arXiv:2607.14450, 16 Jul
2026) and Saad & Ting (arXiv:2603.11015, hierarchical Bayesian reanalysis of the
same 36 systems, γ ≈ 1.12 vs 1.60). Summary reframed as symmetric: the
Newtonian-favouring analyses are disputed on sample construction, the
anomaly-favouring ones on orbital modelling.

**Updated — §12.3 cosmological dependence.** Documents that Θ changed from
2.7×10²⁶ m (2007, 2012) to 8.8×10²⁶ m (2017 onward), moving the predicted floor
from 6.7×10⁻¹⁰ to ~2×10⁻¹⁰ m s⁻² — a factor of ~3.3 in a programme claiming no
adjustable parameters. Also notes the "2" in `2c²/Θ` is `βπ²` with β rounded to
0.2; exact `π²b_W = 1.98779`.

**Updated — §11.2 factorial test.** Added a condition-number ceiling `κ(X) ≤ κmax`
with variance-inflation reporting; the Fowler–Nordheim coupling
`I = I(V, d, φ_wf, β_fe, A_emit)` that makes `∂F/∂V|_I ≈ 0` near-untestable in a
pure field-emission geometry; and a mandatory decoupling stage (independent
gap/bias control, thermionic or photoemission source, or injected beam).

**Mechanical.** Resolved the `η` symbol collision (carrier-transfer map is now
`Φ`); added Θ, η_QI, Φ, M^I_ij, κ(X), |A|/a/a′ and Reading E/R to the notation
dictionary; captioned the Appendix B protocol table and added `\endhead`, closing
the Table 2→4 numbering gap; added the implicit-closure algebra and the numerical
value of `a₀` to Appendix A; six new claim-status ledger rows; new **Criterion 2
(Closure-declaration failure)**.

### `tcgs_bridge_companion.tex` — TCGS bridge

**Added — §6.2, from nonlinear response to observed kernel.** The protocol
treated `K^obs_χ` as a linear operator while the motivating physics is nonlinear.
Now defined as a Fréchet derivative at a declared base point, with a declared
perturbation amplitude, probe design, measured linearity residual `Λ_χ ≤ τ_χ`, and
**Remark 1** limiting the protocol's scope to linearized responses.

**Added — Lemma 1 (shared-sector absorption), with proof.** The claim that the
two shared-sector reconstructions agree under exact factorization is not
automatic: it requires the intersection projectors to be *orthogonal* with
respect to the anchored metric. Self-adjointness is exactly where the metric
enters, which is why it may not be chosen after residuals are inspected.

**Added — Remark 2 (gauge-counting convention).** `D_gauge` was at risk of double
subtraction: `D_r^nom` already quotients the GL(r) factorization redundancy, and
`D^eff_br` is by construction gauge-invariant. Rule stated: `D_gauge = 0` for
frozen maps; otherwise it is the rank of the infinitesimal action on the *joint*
parameter. Notes that the naive bound `N_in² + N_out² = 18` exceeds
`D_r^nom = 8` in the paper's own synthetic design.

**Added — `Σ^bridge_χ` term.** Propagated uncertainty of the fitted bridge
operator, `J_χ Cov(vec K̂_br) J_χᵀ`. Its omission made the residual test
anti-conservative.

**Fixed — symbol collision.** `N_out` was doing duty as both the dimension of the
output bridge space and the intersection of output kernels. The latter is now
`𝒦_out`, named the *portfolio-blind subspace*; both are in the notation dictionary.

**Fixed — counterexample computed on the identifiable image.** §10.3 quoted
`D_r^nom = 196`, contradicting the paper's own definition of `ρ_cap` via
`D^eff_tot`. Now shows the two coincide there (200 conditions on 196 unknowns,
`ker L_P = {0}` generically) and notes that *any* `ρmax < 0.98` rejects the
construction — the margin, not the value 0.25, is the point.

**Fixed — the two nulls separated.** `M_ind` as a *fitted comparison model*
(held-out predictive score) is now distinguished from independent kernels as a
*simulation alternative* (power). Both are required by the protocol; the synthetic
table reports only the second.

**Added.** Estimator disclosure (truncated LS ≠ constrained MLE; size still valid
by Monte Carlo, power is a lower bound); a concrete three-carrier portfolio for
§15.4 (torsion pendulum, cryogenic optomechanical resonator, cold-atom
interferometer, shared low-frequency inertial drive, fourth carrier held back);
five new ledger rows; four new preregistration checklist rows.

**Updated — §12.** The inherited two-closure speculation is withdrawn and
replaced: an unspecified summation domain leaves `𝒜_χ` unfixed by the channel's
physics, so the map can be tuned post hoc through the choice of aperture. That is
the paper's own forbidden item 4 arriving through the input space.

### Bibliography

Added `ChaeYoon2026`, `SaadTing2026`, `McCulloch2011Tajmar`. Restored the missing
issue number on `Banik2024`. Added eprint fields to `ChaeEtAl2026` and
`McCulloch2012`. Flagged the self-cited TCGS preprints for persistent identifiers
(BibTeX has no in-entry comments, so the scaffold sits outside the entry).

### Figures

`scripts/qi_scalar_response.py` regenerates Figure 1 as two panels — Reading E and
Reading R — since the caption now refers to both. The right panel shades the
dynamically inaccessible region and annotates the `F → 0` limit.

---

## v4.0.0 — 22 August 2026

Initial paired release.
