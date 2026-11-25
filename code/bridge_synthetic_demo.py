#!/usr/bin/env python3
"""Reproduce the NEW v6 common-bridge simulation and portfolio-rank checks.

This fully specified experiment does not reconstruct the earlier unavailable
simulation. No seed or parameter is fitted to any previous table.

Run: python code/bridge_synthetic_demo.py [--output-root /path/to/package]
Dependencies: Python >= 3.10, NumPy, Matplotlib. JSON records exact versions,
configuration, maps, truth, individual statistics and binomial uncertainty.

Inputs are iid N(0,1); output noise is iid N(0,0.03^2). Under the null,
Y_c = X_c (P_out,c K P_in,c)^T + noise. Under the alternative, one independent
Gaussian 6x6 kernel is drawn per channel and replicate and held fixed between
calibration and test. Its expected squared Frobenius norm equals the null's.

The declared estimator is common OLS followed by rank-two SVD truncation,
not in general the rank-constrained MLE. Threshold calibration runs through
this exact pipeline, including fit uncertainty, without a chi-square shortcut.

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
    "experiment_id": "bridge_synthetic_v6_new_run_2026-11-25",
    "status": "new_fully_specified_simulation_not_reproduction_of_prior_table",
    "channels": 3,
    "input_dimension_per_channel": 6,
    "output_dimension_per_channel": 6,
    "bridge_input_dimension": 3,
    "bridge_output_dimension": 3,
    "bridge_rank_bound": 2,
    "calibration_inputs_per_channel": 30,
    "heldout_inputs_per_channel": 20,
    "output_noise_standard_deviation": 0.03,
    "inputs": "independent N(0,1) in every coordinate",
    "noise": "independent N(0,sigma^2) in every output coordinate",
    "truth_singular_values": [1.0, 0.5, 0.0],
    "alternative_kernel_entry_variance": 1.25 / 36.0,
    "alternative_expected_squared_frobenius_norm_per_channel": 1.25,
    "null_calibration_replicates": 2500,
    "null_evaluation_replicates": 2500,
    "alternative_evaluation_replicates": 2500,
    "nominal_significance": 0.05,
    "quantile_probability": 0.95,
    "quantile_method": "linear",
    "rejection_rule": "T > calibrated_threshold",
    "random_generator": "numpy.random.Generator(numpy.random.PCG64(seed))",
    "seeds": {
        "channel_maps": 20260822,
        "common_bridge_truth": 20260823,
        "threshold_calibration": 20260824,
        "null_size_evaluation": 20260825,
        "independent_kernel_power_evaluation": 20260826,
        "gaussian_portfolio_rank_checks": 20260827,
    },
    "map_generation_order": "channel index increasing; output map then input map",
    "replicate_draw_order": "alternative kernels if applicable; X_cal; noise_cal; X_test; noise_test",
    "estimator": "common OLS K followed by rank-2 truncated SVD",
    "statistic": "sum over held-out output residuals of residual^2 / sigma^2",
    "interval": "two-sided Wilson 95%, conditional on the frozen MC threshold",
    "wilson_z": 1.959963984540054,
    "demonstration_capacity_ceiling": 0.25,
}


def rng_from_seed(seed):
    return np.random.Generator(np.random.PCG64(seed))


def signed_qr_columns(matrix):
    """Reduced QR with nonnegative diagonal R fixes the sign ambiguity."""
    q, r = np.linalg.qr(matrix, mode="reduced")
    signs = np.where(np.diag(r) < 0.0, -1.0, 1.0)
    return q * signs[None, :]


def rotation(rng, size):
    q = signed_qr_columns(rng.standard_normal((size, size)))
    if np.linalg.det(q) < 0.0:
        q[:, -1] *= -1.0
    return q


def fixed_design():
    rng = rng_from_seed(CONFIG["seeds"]["channel_maps"])
    p_out, p_in = [], []
    for _ in range(3):
        p_out.append(signed_qr_columns(rng.standard_normal((6, 3))))
        p_in.append(signed_qr_columns(rng.standard_normal((6, 3))).T)
    p_out, p_in = np.stack(p_out), np.stack(p_in)
    rng = rng_from_seed(CONFIG["seeds"]["common_bridge_truth"])
    left, right = rotation(rng, 3), rotation(rng, 3)
    k_true = left @ np.diag(CONFIG["truth_singular_values"]) @ right.T
    return p_out, p_in, k_true, p_out @ k_true @ p_in


def fit_common_bridge(x_cal, y_cal, p_out, p_in):
    """K_OLS=solve(sum Z^T Z, sum Z^T Y P_out)^T, then Frobenius SVD truncation."""
    z_cal = x_cal @ np.swapaxes(p_in, -1, -2)
    projected_y = y_cal @ p_out
    gram = np.einsum("cni,cnj->ij", z_cal, z_cal)
    cross = np.einsum("cni,cnj->ij", z_cal, projected_y)
    k_ols = np.linalg.solve(gram, cross).T
    left, singular_values, right_t = np.linalg.svd(k_ols, full_matrices=False)
    r = CONFIG["bridge_rank_bound"]
    return (left[:, :r] * singular_values[:r]) @ right_t[:r]


def simulate(seed, replicates, p_out, p_in, common_kernels, independent_alternative=False):
    """Generate fresh independent calibration and test data in every replicate."""
    rng = rng_from_seed(seed)
    sigma = CONFIG["output_noise_standard_deviation"]
    alt_sd = np.sqrt(CONFIG["alternative_kernel_entry_variance"])
    statistics = np.empty(replicates)
    kernel_norms = np.empty((replicates, 3)) if independent_alternative else None
    for index in range(replicates):
        if independent_alternative:
            kernels = alt_sd * rng.standard_normal((3, 6, 6))
            kernel_norms[index] = np.sum(kernels**2, axis=(1, 2))
        else:
            kernels = common_kernels
        x_cal = rng.standard_normal((3, 30, 6))
        noise_cal = sigma * rng.standard_normal((3, 30, 6))
        x_test = rng.standard_normal((3, 20, 6))
        noise_test = sigma * rng.standard_normal((3, 20, 6))
        kernels_t = np.swapaxes(kernels, -1, -2)
        y_cal, y_test = x_cal @ kernels_t + noise_cal, x_test @ kernels_t + noise_test
        k_fit = fit_common_bridge(x_cal, y_cal, p_out, p_in)
        fitted_kernels = p_out @ k_fit @ p_in
        residual = y_test - x_test @ np.swapaxes(fitted_kernels, -1, -2)
        statistics[index] = np.sum((residual / sigma) ** 2)
    return statistics, kernel_norms


def wilson_interval(successes, trials):
    z = CONFIG["wilson_z"]
    rate = successes / trials
    denominator = 1.0 + z * z / trials
    center = (rate + z * z / (2.0 * trials)) / denominator
    half_width = z * np.sqrt(rate * (1.0 - rate) / trials + z * z / (4.0 * trials**2)) / denominator
    lower = 0.0 if successes == 0 else max(0.0, float(center - half_width))
    upper = 1.0 if successes == trials else min(1.0, float(center + half_width))
    return [lower, upper]


def distribution_summary(values, threshold, independent_evaluation=True):
    rejected = int(np.count_nonzero(values > threshold))
    probs = [0.0, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 1.0]
    quantiles = np.quantile(values, probs, method="linear")
    return {
        "replicates": int(values.size),
        "rejections": rejected,
        "acceptances": int(values.size - rejected),
        "rejection_rate": rejected / values.size,
        "rejection_rate_wilson_95": wilson_interval(rejected, values.size) if independent_evaluation else None,
        "interval_applicability": ("independent evaluation conditional on frozen threshold" if independent_evaluation
                                   else "not applicable: these same replicates determine the quantile threshold"),
        "mean": float(np.mean(values)),
        "standard_deviation_ddof_1": float(np.std(values, ddof=1)),
        "quantiles": {str(p): float(v) for p, v in zip(probs, quantiles)},
        "median_ratio_to_threshold": float(np.median(values) / threshold),
        "minimum_log10_ratio": float(np.min(np.log10(values / threshold))),
        "maximum_log10_ratio": float(np.max(np.log10(values / threshold))),
    }


def numerical_rank_details(matrix):
    singular_values = np.linalg.svd(matrix, compute_uv=False)
    tolerance = float(max(matrix.shape) * np.finfo(float).eps * singular_values[0])
    rank = int(np.count_nonzero(singular_values > tolerance))
    return {
        "matrix_shape": list(matrix.shape),
        "numerical_rank": rank,
        "rank_tolerance": tolerance,
        "smallest_retained_singular_value": float(singular_values[rank - 1]) if rank else None,
        "largest_discarded_singular_value": float(singular_values[rank]) if rank < singular_values.size else None,
        "singular_values": singular_values.tolist(),
    }


def verify_two_channel_portfolios():
    """Two 10x10 observation channels on NxN bridges, N=14 and 18.

    Exact coordinate maps select overlapping sets of 10 coordinates. Rows of
    the portfolio matrix are coordinate basis vectors. Counting distinct rows
    proves rank=2*100-(20-N)^2 exactly, independently of numerical tolerances.
    Separately, iid Gaussian maps verify the generic ranks numerically.
    """
    rng = rng_from_seed(CONFIG["seeds"]["gaussian_portfolio_rank_checks"])
    results = {}
    for n in (14, 18):
        d = 10
        sets = [list(range(d)), list(range(n-d, n))]
        observed_entries, blocks = set(), []
        eye = np.eye(n)
        for indices in sets:
            po, pi = eye[indices, :], eye[:, indices]
            blocks.append(np.kron(pi.T, po))
            observed_entries.update((row, col) for row in indices for col in indices)
        overlap = len(set(sets[0]).intersection(sets[1]))
        exact_rank = len(observed_entries)
        assert exact_rank == 200 - overlap**2 == {14: 164, 18: 196}[n]
        coordinate_details = numerical_rank_details(np.vstack(blocks))
        assert coordinate_details["numerical_rank"] == exact_rank
        gaussian_blocks, gaussian_maps = [], []
        for _ in range(2):
            po, pi = rng.standard_normal((d, n)), rng.standard_normal((n, d))
            gaussian_blocks.append(np.kron(pi.T, po))
            gaussian_maps.append((po, pi))
        gaussian_matrix = np.vstack(gaussian_blocks)
        gaussian_details = numerical_rank_details(gaussian_matrix)
        assert gaussian_details["numerical_rank"] == exact_rank
        probe = rng.standard_normal((n, n))
        direct = np.concatenate([(po @ probe @ pi).reshape(-1, order="F") for po, pi in gaussian_maps])
        kronecker = gaussian_matrix @ probe.reshape(-1, order="F")
        vectorization_error = float(np.max(np.abs(direct - kronecker)))
        assert np.allclose(direct, kronecker, atol=1e-10, rtol=1e-12)
        results[str(n)] = {
            "bridge_input_and_output_dimension": n,
            "channels": 2,
            "channel_input_and_output_dimension": d,
            "rank_bound": n,
            "rank_constraint_is_vacuous": True,
            "nominal_bridge_dimension": n*n,
            "unrestricted_channel_kernel_dimension": 200,
            "generic_input_subspace_intersection_dimension": overlap,
            "generic_output_dual_subspace_intersection_dimension": overlap,
            "coordinate_index_sets_zero_based": sets,
            "exact_coordinate_support_rank": exact_rank,
            "rank_formula": "2*d^2 - max(0,2*d-N)^2",
            "effective_capacity_ratio": exact_rank/200,
            "coordinate_numerical_check": coordinate_details,
            "gaussian_numerical_check": gaussian_details,
            "kronecker_column_vectorization_max_abs_error": vectorization_error,
        }
    return results


def validate_design(p_out, p_in, k_true, kernels):
    identity = np.eye(3)
    out_error = float(np.max(np.abs(np.swapaxes(p_out, -1, -2) @ p_out - identity)))
    in_error = float(np.max(np.abs(p_in @ np.swapaxes(p_in, -1, -2) - identity)))
    assert out_error < 1e-12 and in_error < 1e-12
    # Deterministic probes leave all Monte Carlo streams untouched.
    x = np.arange(3*30*6, dtype=float).reshape(3, 30, 6)
    x = np.sin(x * np.sqrt(2.0)) + np.cos(x * np.sqrt(3.0))
    y = x @ np.swapaxes(kernels, -1, -2)
    recovery_error = float(np.linalg.norm(fit_common_bridge(x, y, p_out, p_in) - k_true))
    assert recovery_error < 1e-10
    # Independently verify the compact normal equations against full vectorized
    # least squares on a nonzero-noise deterministic example.
    perturbation = 0.03 * np.sin(np.arange(y.size).reshape(y.shape) * np.sqrt(5.0))
    y_noisy = y + perturbation
    z = x @ np.swapaxes(p_in, -1, -2)
    full_design = np.vstack([np.kron(row[None, :], p_out[c])
                            for c in range(3) for row in z[c]])
    vectorized_ols = np.linalg.lstsq(full_design, y_noisy.reshape(-1), rcond=None)[0].reshape(3, 3, order="F")
    u, s, vt = np.linalg.svd(vectorized_ols, full_matrices=False)
    reference_fit = (u[:, :2] * s[:2]) @ vt[:2]
    estimator_error = float(np.linalg.norm(reference_fit - fit_common_bridge(x, y_noisy, p_out, p_in)))
    assert estimator_error < 1e-10
    portfolio = np.vstack([np.kron(pi.T, po) for po, pi in zip(p_out, p_in)])
    linear_rank = numerical_rank_details(portfolio)["numerical_rank"]
    assert linear_rank == 9
    return {
        "output_map_orthogonality_max_abs_error": out_error,
        "input_map_orthogonality_max_abs_error": in_error,
        "noiseless_bridge_recovery_frobenius_error": recovery_error,
        "normal_equations_vs_vectorized_least_squares_frobenius_error": estimator_error,
        "fixed_portfolio_linear_rank_on_all_3x3_matrices": linear_rank,
        "identifiable_rank_2_variety_dimension": 8,
        "independent_kernel_capacity": 27,
        "effective_capacity_ratio": 8/27,
        "capacity_denominator": "independent mean kernels conditional on the frozen readable ranges",
        "raw_physical_kernel_capacity": 108,
        "raw_physical_capacity_ratio_diagnostic": 8/108,
        "capacity_gate_passes": 8/27 <= CONFIG["demonstration_capacity_ceiling"],
        "calibration_scalar_observations": 540,
        "informative_projected_calibration_scalar_observations": 270,
        "parameter_to_calibration_data_ratio": 8/270,
        "raw_parameter_to_calibration_data_ratio_diagnostic": 8/540,
        "heldout_scalar_observations": 360,
        "true_bridge_singular_values_computed": np.linalg.svd(k_true, compute_uv=False).tolist(),
        "common_kernel_squared_frobenius_norms": np.sum(kernels**2, axis=(1, 2)).tolist(),
    }


def plot_results(null_values, alternative_values, threshold, target):
    null_log, alt_log = np.log10(null_values/threshold), np.log10(alternative_values/threshold)
    combined = np.concatenate([null_log, alt_log])
    pad = 0.04 * float(np.ptp(combined))
    limits = [float(np.min(combined)-pad), float(np.max(combined)+pad)]
    bins = np.linspace(*limits, 141)
    with plt.rc_context({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42, "ps.fonttype": 42, "axes.labelsize": 11}):
        fig, ax = plt.subplots(figsize=(8.0, 3.75), constrained_layout=True)
        ax.hist(null_log, bins=bins, density=True, color="#236B8E", alpha=0.85,
                label="Common bridge: independent null evaluation", edgecolor="white", linewidth=0.35)
        ax.hist(alt_log, bins=bins, density=True, color="#C65D3B", alpha=0.85,
                label="Independent channel kernels", edgecolor="white", linewidth=0.35)
        ax.axvline(0.0, color="#222222", linewidth=1.15, linestyle="--", label="Rejection threshold")
        ax.set_xlim(limits)
        ax.set_xlabel(r"$\log_{10}(T/q_{0.95})$")
        ax.set_ylabel("Probability density")
        ax.grid(axis="y", color="#D9DDE3", linewidth=0.6, alpha=0.65)
        ax.set_axisbelow(True)
        ax.legend(loc="upper right", frameon=False, fontsize=9)
        ax.set_title("Held-out discrimination in the three-channel control", loc="left", fontsize=12, pad=12)
        fig.savefig(target, format="pdf", metadata={
            "Title": "Held-out discrimination in the three-channel control",
            "Author": "Henry Arellano-Pena",
            "Subject": "All 2500 null and 2500 independent-kernel evaluation statistics; no tail clipping",
            "CreationDate": datetime(2026, 11, 25, tzinfo=timezone.utc),
            "ModDate": datetime(2026, 11, 25, tzinfo=timezone.utc),
        })
        plt.close(fig)
    return {"x_limits": limits, "shared_histogram_bins": 140, "tail_clipping": False,
            "plotted_null_replicates": int(null_values.size), "plotted_alternative_replicates": int(alternative_values.size)}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output-root", type=Path, default=Path(__file__).resolve().parents[1])
    root = parser.parse_args().output_root.resolve()
    result_dir, figure_dir = root/"results", root/"figures"
    result_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    po, pi, k_true, common_kernels = fixed_design()
    verification = validate_design(po, pi, k_true, common_kernels)
    calibration, _ = simulate(CONFIG["seeds"]["threshold_calibration"], 2500, po, pi, common_kernels)
    threshold = float(np.quantile(calibration, 0.95, method="linear"))
    null, _ = simulate(CONFIG["seeds"]["null_size_evaluation"], 2500, po, pi, common_kernels)
    alternative, alternative_norms = simulate(CONFIG["seeds"]["independent_kernel_power_evaluation"],
                                            2500, po, pi, common_kernels, True)
    ranks = verify_two_channel_portfolios()
    figure_path = figure_dir/"bridge_synthetic_residuals.pdf"
    figure_config = plot_results(null, alternative, threshold, figure_path)
    results = {
        "configuration": CONFIG,
        "software": {"python": platform.python_version(), "numpy": np.__version__,
                     "matplotlib": matplotlib.__version__, "platform": platform.platform(),
                     "python_implementation": platform.python_implementation()},
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "design": {"P_out": po.tolist(), "P_in": pi.tolist(), "K_true": k_true.tolist(),
                   "common_channel_kernels": common_kernels.tolist()},
        "design_verification": verification,
        "threshold": threshold,
        "threshold_calibration": distribution_summary(calibration, threshold, independent_evaluation=False),
        "null_size_evaluation": distribution_summary(null, threshold),
        "independent_kernel_power_evaluation": distribution_summary(alternative, threshold),
        "alternative_realized_mean_squared_kernel_norm_per_channel": np.mean(alternative_norms, axis=0).tolist(),
        "alternative_expected_squared_kernel_norm_per_channel": 1.25,
        "statistic_arrays": {"threshold_null": calibration.tolist(), "evaluation_null": null.tolist(),
                             "evaluation_independent_kernels": alternative.tolist()},
        "two_channel_portfolio_rank_checks": ranks,
        "figure": figure_config,
        "interpretation": [
            "New declared v6 simulation, not reconstruction of an earlier table.",
            "Threshold calibrated under one fixed common bridge and fixed independently generated maps.",
            "Evaluation RNG streams independent of threshold calibration and of one another.",
            "Wilson intervals condition on the frozen MC threshold and exclude threshold calibration uncertainty.",
            "Power refers only to the declared Gaussian independent-kernel ensemble.",
            "2500/2500 rejection is empirical power 1, not proof of population power exactly 1.",
            "No power ordering of OLS-plus-SVD relative to constrained least squares is asserted.",
            "V8 capacity accounting conditions on the frozen ranges: 8/27 exceeds the 0.25 ceiling.",
            "The raw ratio 8/108 is retained only as a physical-coordinate diagnostic.",
            "270 informative projected calibration responses exclude known-noise-only output directions; this is an observation count, not a design rank.",
        ],
    }
    result_path = result_dir/"bridge_synthetic_results.json"
    result_path.write_text(json.dumps(results, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    (result_dir/"bridge_portfolio_ranks.json").write_text(json.dumps({
        "seed": CONFIG["seeds"]["gaussian_portfolio_rank_checks"], "checks": ranks},
        indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(json.dumps({"threshold": threshold, "null_size": results["null_size_evaluation"],
        "independent_kernel_power": results["independent_kernel_power_evaluation"],
        "portfolio_ranks": {k: v["exact_coordinate_support_rank"] for k, v in ranks.items()},
        "results": str(result_path), "figure": str(figure_path)}, indent=2))


if __name__ == "__main__":
    main()
