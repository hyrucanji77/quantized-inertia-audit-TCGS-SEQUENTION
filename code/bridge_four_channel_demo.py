#!/usr/bin/env python3
"""V8: a four-channel control with capacity conditional on frozen map ranges.

Run, in this order, from any directory:
    python code/bridge_four_channel_demo.py --design-only
    python code/bridge_four_channel_demo.py

The first command freezes the complete computational design before simulation;
it is not a public preregistration. The second checks the configuration exactly,
checks regenerated arrays with numerical tolerances, and uses the saved arrays.
Computed verification diagnostics are checked against tolerances, not for bit
identity. An implementation-only loader amendment preserves the original design,
result records and generating script; see the separate verification record.
All three Monte Carlo ensembles have 2,500 replicates. No run is selected or
tuned to match another result. Existing three-channel controls are untouched.

Four physical 6x6 kernels are read through fixed rank-three input/output maps.
The conditional independent-model dimension is therefore 4*3*3=36, not 144.
A common rank-two 3x3 bridge has visible dimension eight. The common estimator
is unconstrained pooled OLS followed by Frobenius rank-two SVD truncation, not
in general the rank-constrained MLE. The full estimator/test pipeline is used
in Monte Carlo threshold calibration. No nominal chi-square approximation is
used. The alternative retains every frozen map range and each local rank bound;
its four readable blocks have independent orientations and fixed singular
values (1,.5,0), and remain fixed across calibration and holdout in a replicate.

An independent unrestricted 3x3 readable-block model is fitted on calibration
data for a descriptive predictive-score comparison. That comparison does not
establish the physical protocol's evidence, novelty, or standard-physics gates.
Size/power and Wilson intervals are pointwise for the declared frozen truth and
alternative ensemble, conditional on the realized Monte Carlo threshold.

Copyright (c) 2026 Henry Arellano-Pena. CC BY 4.0.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


CONFIG = {
    "experiment_id": "v8_four_channel_conditional_capacity_2026-11-25",
    "declaration_status": "computational specification saved before Monte Carlo; not public preregistration",
    "channels": 4,
    "physical_input_dimension_per_channel": 6,
    "physical_output_dimension_per_channel": 6,
    "bridge_input_dimension": 3,
    "bridge_output_dimension": 3,
    "bridge_rank_bound": 2,
    "map_rank_per_channel": 3,
    "calibration_inputs_per_channel": 30,
    "heldout_inputs_per_channel": 20,
    "output_noise_standard_deviation": 0.03,
    "inputs": "iid N(0,I6), independent between channels, splits and replicates",
    "noise": "iid N(0,sigma^2) in every physical output coordinate",
    "common_truth_singular_values": [1.0, 0.5, 0.0],
    "alternative_singular_values_per_readable_block": [1.0, 0.5, 0.0],
    "alternative_squared_frobenius_norm_per_readable_block": 1.25,
    "alternative": "independent SO(3) left/right orientations per channel and replicate; readable blocks held fixed across splits",
    "alternative_scope": "independent rank-two readable blocks inside all fixed map ranges; ensemble-specific, not local or worst-case power",
    "threshold_replicates": 2500,
    "null_size_replicates": 2500,
    "alternative_power_replicates": 2500,
    "nominal_significance": 0.05,
    "quantile_probability": 0.95,
    "quantile_method": "linear",
    "rejection_rule": "common_physical_heldout_T > frozen_MC_threshold",
    "common_estimator": "pooled unconstrained OLS K, followed by Frobenius rank-two SVD truncation",
    "common_statistic": "physical held-out sum of squared residuals divided by known sigma^2",
    "conditional_comparator": "independent unrestricted 3x3 readable blocks fitted by calibration OLS; dimension 36",
    "score_difference": "T_independent - T_common; positive means better common held-out Gaussian log score; log-score difference is half this value",
    "comparison_scope": "descriptive only; no evidence gate, standard-physics test or complete-protocol success claimed",
    "rho_max": 0.25,
    "eta_max": 0.1,
    "capacity_denominator": "36: independent kernels conditional on all frozen readable ranges",
    "capacity_interpretation": "combined common-bridge and rank-two restrictions relative to independent unrestricted readable blocks",
    "secondary_capacity_reference": "independent rank-two readable blocks have dimension 4*8=32; 8/32=.25, reported without changing the declared comparator",
    "data_denominator": "360: calibration scalar outputs retained by the rank-three output maps; sigma and maps are known",
    "size_scope": "pointwise at the frozen common truth; not uniform over the composite null",
    "interval": "pointwise two-sided Wilson 95%, conditional on the realized frozen MC threshold; excludes threshold uncertainty",
    "wilson_z": 1.959963984540054,
    "random_generator": "numpy.random.Generator(numpy.random.PCG64(seed))",
    "seeds": {
        "channel_maps": 2026112801,
        "common_bridge_truth": 2026112802,
        "threshold_calibration": 2026112803,
        "null_size_evaluation": 2026112804,
        "alternative_power_evaluation": 2026112805,
        "independent_verification_probes": 2026112806,
    },
    "map_draw_order": "channel index increasing; output 6x3 isometry then input 6x3 isometry, transposed",
    "rotation_rule": "signed QR with nonnegative R diagonal, then flip final column if determinant is negative",
    "alternative_draw_order": "channel index increasing; left SO(3) orientation then right SO(3) orientation",
    "replicate_draw_order": "alternative blocks if applicable; X_cal; noise_cal; X_test; noise_test",
    "figure": "shared 160-bin histograms of log10(T/threshold), all independent evaluation replicates, no tail clipping",
}


def rng(seed):
    return np.random.Generator(np.random.PCG64(seed))


def signed_qr(matrix):
    q, r = np.linalg.qr(matrix, mode="reduced")
    return q * np.where(np.diag(r) < 0.0, -1.0, 1.0)[None, :]


def rotation(generator, n):
    q = signed_qr(generator.standard_normal((n, n)))
    if np.linalg.det(q) < 0.0:
        q[:, -1] *= -1.0
    return q


def construct_design():
    generator = rng(CONFIG["seeds"]["channel_maps"])
    output_maps, input_maps = [], []
    for _ in range(CONFIG["channels"]):
        output_maps.append(signed_qr(generator.standard_normal((6, 3))))
        input_maps.append(signed_qr(generator.standard_normal((6, 3))).T)
    generator = rng(CONFIG["seeds"]["common_bridge_truth"])
    left, right = rotation(generator, 3), rotation(generator, 3)
    truth = left @ np.diag(CONFIG["common_truth_singular_values"]) @ right.T
    output_maps, input_maps = np.stack(output_maps), np.stack(input_maps)
    return {"P_out": output_maps, "P_in": input_maps, "K_true": truth,
            "truth_left_orientation": left, "truth_right_orientation": right,
            "common_physical_kernels": output_maps @ truth @ input_maps}


def rank_details(matrix):
    s = np.linalg.svd(matrix, compute_uv=False)
    tol = float(max(matrix.shape) * np.finfo(float).eps * s[0])
    return {"shape": list(matrix.shape), "rank": int(np.count_nonzero(s > tol)),
            "tolerance": tol, "singular_values": s.tolist()}


def truncate_rank_two(matrix):
    u, s, vt = np.linalg.svd(matrix, full_matrices=False)
    return (u[:, :2] * s[:2]) @ vt[:2]


def fit_models(x, y, design):
    z = x @ design["P_in"].transpose(0, 2, 1)
    projected_y = y @ design["P_out"]
    gram = z.transpose(0, 2, 1) @ z
    cross = z.transpose(0, 2, 1) @ projected_y
    common_ols = np.linalg.solve(np.sum(gram, axis=0), np.sum(cross, axis=0)).T
    common_rank_two = truncate_rank_two(common_ols)
    independent = np.linalg.solve(gram, cross).transpose(0, 2, 1)
    return common_rank_two, independent, common_ols


def draw_alternative(generator):
    return np.stack([rotation(generator, 3) @ np.diag(CONFIG["alternative_singular_values_per_readable_block"])
                     @ rotation(generator, 3).T for _ in range(CONFIG["channels"])])


def validate_design(design):
    po, pi, truth = design["P_out"], design["P_in"], design["K_true"]
    out_error = float(np.max(np.abs(po.transpose(0, 2, 1) @ po - np.eye(3))))
    in_error = float(np.max(np.abs(pi @ pi.transpose(0, 2, 1) - np.eye(3))))
    assert out_error < 1e-12 and in_error < 1e-12
    operator = np.vstack([np.kron(p.T, q) for q, p in zip(po, pi)])
    linear = rank_details(operator)
    assert linear["rank"] == 9
    # Tangent to the rank-two variety: in the truth's singular-vector basis,
    # all entries except the bottom-right coordinate can vary to first order.
    left, right = design["truth_left_orientation"], design["truth_right_orientation"]
    tangent = np.column_stack([np.outer(left[:, i], right[:, j]).reshape(-1, order="F")
                               for i in range(3) for j in range(3) if i < 2 or j < 2])
    tangent_image = rank_details(operator @ tangent)
    assert tangent_image["rank"] == 8
    generator = rng(CONFIG["seeds"]["independent_verification_probes"])
    x = generator.standard_normal((4, 30, 6))
    noiseless_y = x @ design["common_physical_kernels"].transpose(0, 2, 1)
    fit, independent, _ = fit_models(x, noiseless_y, design)
    recovery = float(np.linalg.norm(fit - truth))
    independent_recovery = float(np.max(np.abs(independent - truth)))
    assert recovery < 1e-10 and independent_recovery < 1e-10
    y = noiseless_y + 0.03 * generator.standard_normal(noiseless_y.shape)
    fit, local, pooled_ols = fit_models(x, y, design)
    z = x @ pi.transpose(0, 2, 1)
    full_design = np.vstack([np.kron(row[None, :], po[c]) for c in range(4) for row in z[c]])
    reference = np.linalg.lstsq(full_design, y.reshape(-1), rcond=None)[0].reshape(3, 3, order="F")
    ols_error = float(np.linalg.norm(reference - pooled_ols))
    truncation_error = float(np.linalg.norm(truncate_rank_two(reference) - fit))
    local_errors = []
    for c in range(4):
        a = np.vstack([np.kron(row[None, :], po[c]) for row in z[c]])
        ref = np.linalg.lstsq(a, y[c].reshape(-1), rcond=None)[0].reshape(3, 3, order="F")
        local_errors.append(float(np.linalg.norm(ref - local[c])))
    assert max([ols_error, truncation_error] + local_errors) < 1e-10
    blocks = draw_alternative(generator)
    kernels = po @ blocks @ pi
    recovered = po.transpose(0, 2, 1) @ kernels @ pi.transpose(0, 2, 1)
    range_error = float(np.max(np.abs(recovered - blocks)))
    singular_values = np.linalg.svd(blocks, compute_uv=False)
    assert range_error < 1e-12 and np.max(np.abs(singular_values[:, 2])) < 1e-12
    assert np.max(np.abs(np.sum(blocks**2, axis=(1, 2)) - 1.25)) < 1e-12
    return {
        "output_map_orthogonality_max_abs_error": out_error,
        "input_map_orthogonality_max_abs_error": in_error,
        "full_linear_observation_operator": linear,
        "rank_two_tangent_image": tangent_image,
        "finite_bridge_observation_kernel_dimension": 0,
        "conditional_independent_kernel_capacity": 36,
        "common_bridge_visible_capacity": 8,
        "conditional_capacity_ratio": 8 / 36,
        "secondary_independent_rank_two_readable_capacity": 32,
        "secondary_commonality_only_capacity_ratio": 8 / 32,
        "rho_max": CONFIG["rho_max"],
        "capacity_gate_passes": bool(8 / 36 <= CONFIG["rho_max"]),
        "unrestricted_physical_kernel_capacity_diagnostic_only": 144,
        "raw_capacity_ratio_diagnostic_only": 8 / 144,
        "calibration_raw_scalar_outputs": 720,
        "calibration_informative_scalar_outputs": 360,
        "calibration_count_interpretation": "independent informative scalar output observations conditional on inputs/maps; not Fisher-information or design-matrix rank",
        "parameter_to_informative_calibration_data_ratio": 8 / 360,
        "eta_max": CONFIG["eta_max"],
        "data_gate_passes": bool(8 / 360 <= CONFIG["eta_max"]),
        "heldout_raw_scalar_outputs": 480,
        "heldout_retained_scalar_outputs": 240,
        "noiseless_common_recovery_frobenius_error": recovery,
        "noiseless_independent_recovery_max_abs_error": independent_recovery,
        "pooled_OLS_vs_vectorized_least_squares_frobenius_error": ols_error,
        "truncated_pooled_OLS_vs_vectorized_frobenius_error": truncation_error,
        "independent_OLS_vs_vectorized_frobenius_errors": local_errors,
        "verification_alternative_readable_range_recovery_max_abs_error": range_error,
        "verification_alternative_singular_values": singular_values.tolist(),
        "common_truth_singular_values_computed": np.linalg.svd(truth, compute_uv=False).tolist(),
    }


def jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {key: jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    return value


def frozen_payload(design, verification):
    return {"configuration": CONFIG, "fixed_design": jsonable(design),
            "design_verification": verification,
            "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


def freeze_design(path, payload):
    if path.exists():
        raise FileExistsError("A saved design must be loaded and verified, never overwritten.")
    content = {"design_saved_utc": datetime.now(timezone.utc).isoformat(),
               "frozen_payload": payload}
    path.write_text(json.dumps(content, indent=2) + "\n")


def load_frozen_design(path, regenerated):
    """Keep saved arrays authoritative; do not require QR/SVD bit identity."""
    saved = json.loads(path.read_text())
    payload = saved["frozen_payload"]
    assert payload["configuration"] == CONFIG, "Configuration changed after freezing; do not run."
    assert set(payload["fixed_design"]) == set(regenerated), "Fixed-design keys changed."
    design, differences = {}, {}
    for key, expected in regenerated.items():
        value = np.asarray(payload["fixed_design"][key], dtype=float)
        assert value.shape == expected.shape, f"Saved array shape differs: {key}"
        assert np.all(np.isfinite(value)), f"Saved array contains nonfinite values: {key}"
        assert np.allclose(value, expected, rtol=2e-12, atol=2e-12), f"Saved array differs: {key}"
        design[key] = value
        differences[key] = float(np.max(np.abs(value - expected)))
    return saved, design, {"relative_tolerance": 2e-12, "absolute_tolerance": 2e-12,
                           "saved_arrays_used_for_simulation": True,
                           "maximum_absolute_differences_by_array": differences}


def simulate(seed, count, design, alternative=False):
    generator = rng(seed)
    po, pi = design["P_out"], design["P_in"]
    sigma = CONFIG["output_noise_standard_deviation"]
    values = {key: np.empty(count) for key in (
        "common_physical_T", "independent_physical_T", "independent_minus_common_T",
        "common_retained_T", "independent_retained_T", "orthogonal_output_noise_T")}
    all_blocks = np.empty((count, 4, 3, 3)) if alternative else None
    for j in range(count):
        if alternative:
            blocks = draw_alternative(generator)
            all_blocks[j] = blocks
            kernels = po @ blocks @ pi
        else:
            kernels = design["common_physical_kernels"]
        x_cal = generator.standard_normal((4, 30, 6))
        noise_cal = sigma * generator.standard_normal((4, 30, 6))
        x_test = generator.standard_normal((4, 20, 6))
        noise_test = sigma * generator.standard_normal((4, 20, 6))
        y_cal = x_cal @ kernels.transpose(0, 2, 1) + noise_cal
        y_test = x_test @ kernels.transpose(0, 2, 1) + noise_test
        common, independent, _ = fit_models(x_cal, y_cal, design)
        common_residual = y_test - x_test @ (po @ common @ pi).transpose(0, 2, 1)
        independent_residual = y_test - x_test @ (po @ independent @ pi).transpose(0, 2, 1)
        common_t = float(np.sum((common_residual / sigma)**2))
        independent_t = float(np.sum((independent_residual / sigma)**2))
        common_retained = float(np.sum((common_residual @ po / sigma)**2))
        independent_retained = float(np.sum((independent_residual @ po / sigma)**2))
        orthogonal = y_test - (y_test @ po) @ po.transpose(0, 2, 1)
        orthogonal_t = float(np.sum((orthogonal / sigma)**2))
        assert np.isclose(common_t, common_retained + orthogonal_t, rtol=1e-11, atol=1e-9)
        assert np.isclose(independent_t, independent_retained + orthogonal_t, rtol=1e-11, atol=1e-9)
        values["common_physical_T"][j] = common_t
        values["independent_physical_T"][j] = independent_t
        values["independent_minus_common_T"][j] = independent_t - common_t
        values["common_retained_T"][j] = common_retained
        values["independent_retained_T"][j] = independent_retained
        values["orthogonal_output_noise_T"][j] = orthogonal_t
    return values, all_blocks


def wilson(successes, total):
    z = CONFIG["wilson_z"]
    rate, denom = successes / total, 1.0 + z**2 / total
    center = (rate + z**2 / (2 * total)) / denom
    half = z * np.sqrt(rate * (1 - rate) / total + z**2 / (4 * total**2)) / denom
    return [0.0 if successes == 0 else float(max(0.0, center - half)),
            1.0 if successes == total else float(min(1.0, center + half))]


def summary(values):
    probabilities = [0.0, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 1.0]
    return {"replicates": int(values.size), "mean": float(np.mean(values)),
            "standard_deviation_ddof_1": float(np.std(values, ddof=1)),
            "quantiles": {str(p): float(v) for p, v in zip(probabilities, np.quantile(values, probabilities, method="linear"))}}


def test_summary(values, threshold, independent=True):
    result = summary(values)
    rejected = int(np.count_nonzero(values > threshold))
    result.update({"rejections": rejected, "rejection_rate": rejected / values.size,
                   "wilson_95": wilson(rejected, values.size) if independent else None,
                   "interval_scope": CONFIG["interval"] if independent else "not applicable: threshold-estimation sample"})
    return result


def plot_results(null, alternative, threshold, destination):
    null_x, alt_x = np.log10(null / threshold), np.log10(alternative / threshold)
    all_x = np.concatenate([null_x, alt_x])
    pad = 0.04 * float(np.ptp(all_x))
    bounds = [float(all_x.min() - pad), float(all_x.max() + pad)]
    bins = np.linspace(*bounds, 161)
    with plt.rc_context({"font.size": 10, "pdf.fonttype": 42, "ps.fonttype": 42,
                         "axes.spines.top": False, "axes.spines.right": False}):
        fig, ax = plt.subplots(figsize=(8.0, 3.8), constrained_layout=True)
        ax.hist(null_x, bins=bins, density=True, color="#246B8E", alpha=0.88,
                edgecolor="white", linewidth=0.3, label="Common rank-two bridge")
        ax.hist(alt_x, bins=bins, density=True, color="#B85E3D", alpha=0.88,
                edgecolor="white", linewidth=0.3, label="Independent rank-two readable blocks")
        ax.axvline(0.0, color="#242A30", linestyle="--", linewidth=1.1, label="Frozen 95% rejection threshold")
        ax.set_xlim(bounds)
        ax.set_xlabel(r"$\log_{10}(T/q_{0.95})$")
        ax.set_ylabel("Probability density")
        ax.set_title("Four-channel control within frozen readable ranges", loc="left", fontsize=12, pad=11)
        ax.grid(axis="y", color="#D9DDE3", linewidth=0.6, alpha=0.6)
        ax.set_axisbelow(True)
        ax.legend(frameon=False, fontsize=9, loc="upper right")
        fig.savefig(destination, format="pdf", metadata={
            "Title": "Four-channel control within frozen readable ranges",
            "Author": "Henry Arellano-Pena",
            "Subject": "Declared four-channel simulation; 2500 null and 2500 rank-two in-range alternative evaluations; no tail clipping",
            "CreationDate": datetime(2026, 11, 25, tzinfo=timezone.utc),
            "ModDate": datetime(2026, 11, 25, tzinfo=timezone.utc),
        })
        plt.close(fig)
    return {"x_limits": bounds, "shared_histogram_bins": 160, "tail_clipping": False,
            "null_replicates": int(null.size), "alternative_replicates": int(alternative.size)}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--design-only", action="store_true")
    args = parser.parse_args()
    root = args.output_root.resolve()
    result_dir, figure_dir = root / "results", root / "figures"
    result_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    regenerated, design_path = construct_design(), result_dir / "bridge_four_channel_design.json"
    current_script_sha256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if args.design_only:
        if design_path.exists():
            stored_design, design, _ = load_frozen_design(design_path, regenerated)
            validate_design(design)
        else:
            verification = validate_design(regenerated)
            freeze_design(design_path, frozen_payload(regenerated, verification))
            stored_design = json.loads(design_path.read_text())
        print(json.dumps({"design_saved": str(design_path), "script_sha256": current_script_sha256,
                          "frozen_script_sha256": stored_design["frozen_payload"]["script_sha256"],
                          "rho_conditional": 8/36, "eta_informative": 8/360}))
        return
    if not design_path.exists():
        parser.error("Run --design-only first, before any Monte Carlo simulation.")
    stored_design, design, regeneration_check = load_frozen_design(design_path, regenerated)
    verification = validate_design(design)
    started_utc = datetime.now(timezone.utc).isoformat()
    calibration, _ = simulate(CONFIG["seeds"]["threshold_calibration"], CONFIG["threshold_replicates"], design)
    threshold = float(np.quantile(calibration["common_physical_T"], CONFIG["quantile_probability"], method=CONFIG["quantile_method"]))
    null, _ = simulate(CONFIG["seeds"]["null_size_evaluation"], CONFIG["null_size_replicates"], design)
    alternative, blocks = simulate(CONFIG["seeds"]["alternative_power_evaluation"], CONFIG["alternative_power_replicates"], design, True)
    alternative_s = np.linalg.svd(blocks, compute_uv=False)
    alternative_norm = np.sum(blocks**2, axis=(2, 3))
    assert np.max(np.abs(alternative_s - [1.0, 0.5, 0.0])) < 1e-12
    assert np.max(np.abs(alternative_norm - 1.25)) < 1e-12
    figure_path = figure_dir / "bridge_four_channel_results.pdf"
    figure = plot_results(null["common_physical_T"], alternative["common_physical_T"], threshold, figure_path)
    result = {
        "configuration": CONFIG, "design_saved_utc": stored_design["design_saved_utc"],
        "simulation_started_utc": started_utc, "simulation_completed_utc": datetime.now(timezone.utc).isoformat(),
        "fixed_design": jsonable(design), "design_verification": verification,
        "script_sha256": current_script_sha256,
        "frozen_script_sha256": stored_design["frozen_payload"]["script_sha256"],
        "fixed_design_regeneration_check": regeneration_check,
        "frozen_design_sha256": hashlib.sha256(design_path.read_bytes()).hexdigest(),
        "software": {"python": platform.python_version(), "python_implementation": platform.python_implementation(),
                     "numpy": np.__version__, "matplotlib": matplotlib.__version__, "platform": platform.platform()},
        "threshold": threshold,
        "threshold_sample": test_summary(calibration["common_physical_T"], threshold, False),
        "null_size": test_summary(null["common_physical_T"], threshold),
        "alternative_power": test_summary(alternative["common_physical_T"], threshold),
        "descriptive_predictive_comparison": {
            "scope": CONFIG["comparison_scope"], "sign_convention": CONFIG["score_difference"],
            "null_independent_minus_common_T": summary(null["independent_minus_common_T"]),
            "alternative_independent_minus_common_T": summary(alternative["independent_minus_common_T"]),
            "null_independent_T": summary(null["independent_physical_T"]),
            "alternative_independent_T": summary(alternative["independent_physical_T"]),
        },
        "alternative_realization_verification": {
            "maximum_singular_value_deviation_from_declared_values": float(np.max(np.abs(alternative_s - [1.0, 0.5, 0.0]))),
            "maximum_squared_norm_deviation_from_1p25": float(np.max(np.abs(alternative_norm - 1.25))),
            "local_ranks_two_by_construction": True, "all_blocks_within_frozen_map_ranges_by_construction": True,
        },
        "per_replicate_statistics": {"threshold_sample": jsonable(calibration), "null_evaluation": jsonable(null),
                                     "alternative_evaluation": jsonable(alternative)},
        "alternative_readable_blocks_per_replicate": blocks.tolist(),
        "figure": {"path": str(figure_path.relative_to(root)), **figure},
    }
    result_path = result_dir / "bridge_four_channel_results.json"
    result_path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: result[key] for key in ("threshold", "null_size", "alternative_power", "descriptive_predictive_comparison")}, indent=2))


if __name__ == "__main__":
    main()
