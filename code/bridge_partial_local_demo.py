#!/usr/bin/env python3
"""V7: partial bridge readout, correlated shared-sector tests, local alternatives.

Run from any directory:
    python code/bridge_partial_local_demo.py --design-only
    python code/bridge_partial_local_demo.py

Only NumPy and Matplotlib are required. All statistical choices below were
specified before the Monte Carlo run. Saving a computational specification is
not public preregistration. This program does not change the v6 extreme control.

Each channel reads only three dimensions of a six-dimensional bridge. The
sectors share two dimensions; the sixth coordinate is wholly unread. OLS on
20 held-out probes estimates each readable block. Partial Moore-Penrose
reconstruction and shared-sector projection yield three correlated pairwise
contrasts. Their joint covariance has rank eight. The global statistic uses
the equivalent eight independent coordinates with their complete covariance.
It tests common compatibility, not the entire global rank-two hypothesis.

Local truth blocks remain rank two and inside the fixed channel-map ranges.
Delta is exact relative RMS response displacement from the fixed baseline,
not distance to the entire compatible family. Power pertains to the declared
random perturbation ensemble; it is not a worst-case separation guarantee.

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
from matplotlib.colors import ListedColormap, BoundaryNorm
import numpy as np


DELTAS = [0.0, 0.005, 0.01, 0.02, 0.027, 0.03, 0.045, 0.063, 0.08, 0.1, 0.15, 0.2]
CONFIG = {
    "experiment_id": "v7_partial_shared_sector_local_power_2026-11-25",
    "declaration_status": "computational specification saved before Monte Carlo; not public preregistration",
    "bridge_input_dimension": 6,
    "bridge_output_dimension": 6,
    "channels": 3,
    "input_dimension_per_channel": 6,
    "output_dimension_per_channel": 6,
    "map_rank_per_channel": 3,
    "readable_coordinate_sets_zero_based": [[0, 1, 2], [0, 1, 3], [0, 1, 4]],
    "shared_input_and_output_dimension": 2,
    "completely_unread_coordinate_zero_based": 5,
    "common_truth_rank": 2,
    "common_truth_singular_values": [1.0, 0.5, 0.0, 0.0, 0.0, 0.0],
    "calibration_samples_per_channel": 30,
    "heldout_samples_per_channel": 20,
    "output_noise_standard_deviation": 0.03,
    "input_law": "iid N(0,1) coordinates, independent between channels and splits",
    "noise_law": "iid N(0,0.03^2) coordinates, independent between channels and splits",
    "null_threshold_replicates": 2500,
    "null_size_replicates": 2500,
    "local_replicates_per_delta": 2500,
    "deltas": DELTAS,
    "significance": 0.05,
    "quantile_probability": 0.95,
    "quantile_method": "linear",
    "rejection_rule": "global_shared_T > frozen_MC_threshold",
    "power_planning_target": 0.95,
    "planning_rule": "first tested positive delta with empirical rejection rate >= 0.95; no interpolation",
    "delta_definition": "sqrt(sum_c ||B_c(delta)-B_c(0)||_F^2 / sum_c ||B_c(0)||_F^2)",
    "alternative": "H_c(delta)=H_c0+delta*lambda*U_c*D_c; D_c iid N(0,1) in R^(2x3)",
    "normalization": "lambda=sqrt(sum ||H_c0||_F^2 / sum ||U_c D_c||_F^2)",
    "alternative_basis": "first two left singular vectors of H_c0; each column sign fixed by its largest absolute entry",
    "alternative_scope": "rank-2 local blocks, fixed channel ranges, ensemble-specific; delta is not distance to entire null family",
    "alternative_draw_frequency": "fresh directions each replicate; kernels held fixed between calibration and heldout splits",
    "shared_test_fit": "unconstrained OLS of each 3x3 readable block on heldout probes, with exact conditional covariance",
    "shared_test_scope": "common-compatible shared sector only; no full global rank-2 adequacy or physical ECL test",
    "pair_order_zero_based": [[0, 1], [0, 2], [1, 2]],
    "global_contrast": "u=(vec(Hshared_0-Hshared_2), vec(Hshared_1-Hshared_2)), column vectorization",
    "global_covariance": "W=[[V0+V2,V2],[V2,V1+V2]], Vc=sigma^2*kron((Zc^T Zc)^(-1)[0:2,0:2],I2)",
    "calibration_predictive_diagnostic": "19-parameter common-compatible quotient vs 27-parameter independent readable blocks; descriptive only",
    "wilson_z": 1.959963984540054,
    "interval_scope": "pointwise 95% Wilson intervals conditional on frozen MC threshold; not simultaneous or threshold-uncertainty intervals",
    "rng": "numpy.random.Generator(numpy.random.PCG64(seed))",
    "seeds": {
        "maps": 2026112501,
        "common_truth": 2026112502,
        "threshold": 2026112503,
        "size": 2026112504,
        "deterministic_verification_probes": 2026112505,
        "local_by_delta": {str(delta): 2026112600+i for i, delta in enumerate(DELTAS)},
    },
    "map_draw_order": "ambient output basis, ambient input basis; then output and input 6x3 isometries for each channel",
    "replicate_draw_order": "local D if applicable; X_cal; noise_cal; X_test; noise_test",
}


def rng(seed):
    return np.random.Generator(np.random.PCG64(seed))


def signed_qr(matrix):
    q, r = np.linalg.qr(matrix, mode="reduced")
    return q * np.where(np.diag(r) < 0, -1.0, 1.0)[None, :]


def rotation(generator, n):
    q = signed_qr(generator.standard_normal((n, n)))
    if np.linalg.det(q) < 0:
        q[:, -1] *= -1
    return q


def rank_info(matrix):
    s = np.linalg.svd(matrix, compute_uv=False)
    tolerance = max(matrix.shape) * np.finfo(float).eps * s[0]
    rank = int(np.count_nonzero(s > tolerance))
    return {"shape": list(matrix.shape), "rank": rank, "tolerance": float(tolerance),
            "singular_values": s.tolist()}


def construct_design():
    generator = rng(CONFIG["seeds"]["maps"])
    q_out, q_in = rotation(generator, 6), rotation(generator, 6)
    local_out, local_in, e_out, e_in = [], [], [], []
    sectors = CONFIG["readable_coordinate_sets_zero_based"]
    for coordinates in sectors:
        local_out.append(signed_qr(generator.standard_normal((6, 3))))
        local_in.append(signed_qr(generator.standard_normal((6, 3))))
        e_out.append(q_out[:, coordinates])
        e_in.append(q_in[:, coordinates])
    local_out, local_in = np.stack(local_out), np.stack(local_in)
    e_out, e_in = np.stack(e_out), np.stack(e_in)
    p_out = local_out @ e_out.transpose(0, 2, 1)
    p_in = e_in @ local_in.transpose(0, 2, 1)
    truth_rng = rng(CONFIG["seeds"]["common_truth"])
    left, right = rotation(truth_rng, 6), rotation(truth_rng, 6)
    factor_a, factor_b = left[:, :2] * np.sqrt([1.0, 0.5]), right[:, :2] * np.sqrt([1.0, 0.5])
    k_true = factor_a @ factor_b.T
    h0 = e_out.transpose(0, 2, 1) @ k_true @ e_in
    perturbation_bases = []
    for h in h0:
        u, _, _ = np.linalg.svd(h, full_matrices=False)
        u = u[:, :2]
        for col in range(2):
            pivot = np.argmax(np.abs(u[:, col]))
            if u[pivot, col] < 0:
                u[:, col] *= -1
        perturbation_bases.append(u)
    visible_entries = sorted({(row, col) for sector in sectors for row in sector for col in sector})
    entry_index = {entry: index for index, entry in enumerate(visible_entries)}
    block_maps = np.zeros((3, 9, len(visible_entries)))
    for c, sector in enumerate(sectors):
        for row in range(3):
            for col in range(3):
                block_maps[c, row+3*col, entry_index[(sector[row], sector[col])]] = 1
    return {"Q_out": q_out, "Q_in": q_in, "local_out": local_out, "local_in": local_in,
            "E_out": e_out, "E_in": e_in, "P_out": p_out, "P_in": p_in,
            "common_out": q_out[:, :2], "common_in": q_in[:, :2],
            "factor_A": factor_a, "factor_B": factor_b, "K_true": k_true,
            "H0": h0, "perturbation_bases": np.stack(perturbation_bases),
            "visible_entries": visible_entries, "block_maps": block_maps}


def fit_readable_blocks(x, y, design):
    z = x @ design["local_in"]
    projected_y = y @ design["local_out"]
    gram = z.transpose(0, 2, 1) @ z
    cross = z.transpose(0, 2, 1) @ projected_y
    h = np.linalg.solve(gram, cross).transpose(0, 2, 1)
    gram_inv = np.linalg.solve(gram, np.broadcast_to(np.eye(3), gram.shape))
    return h, gram_inv, z, projected_y, gram


def pair_contrast_matrix():
    identity, zero = np.eye(4), np.zeros((4, 4))
    return np.block([[identity, -identity, zero], [identity, zero, -identity], [zero, identity, -identity]])


def shared_statistics(h, gram_inv, design, full_check=False):
    # Explicit partial reconstructions, rather than treating any channel as an inverse.
    reconstructed = design["E_out"] @ h @ design["E_in"].transpose(0, 2, 1)
    shared = design["common_out"].T @ reconstructed @ design["common_in"]
    hs = np.stack([block.reshape(4, order="F") for block in shared])
    sigma2 = CONFIG["output_noise_standard_deviation"]**2
    covariances = np.stack([sigma2*np.kron(g[:2, :2], np.eye(2)) for g in gram_inv])
    v0, v1, v2 = covariances
    independent = np.concatenate([hs[0]-hs[2], hs[1]-hs[2]])
    covariance = np.block([[v0+v2, v2], [v2, v1+v2]])
    global_t = float(independent @ np.linalg.solve(covariance, independent))
    pair_residuals = [hs[a]-hs[b] for a, b in ((0, 1), (0, 2), (1, 2))]
    pair_t = [float(d @ np.linalg.solve(covariances[a]+covariances[b], d))
              for d, (a, b) in zip(pair_residuals, ((0, 1), (0, 2), (1, 2)))]
    details = None
    if full_check:
        block_cov = np.zeros((12, 12))
        for c in range(3):
            block_cov[c*4:(c+1)*4, c*4:(c+1)*4] = covariances[c]
        contrast = pair_contrast_matrix()
        full_cov = contrast @ block_cov @ contrast.T
        full_residual = np.concatenate(pair_residuals)
        redundant_t = float(full_residual @ np.linalg.pinv(full_cov, hermitian=True) @ full_residual)
        assert np.isclose(global_t, redundant_t, atol=1e-9, rtol=1e-10)
        assert rank_info(full_cov)["rank"] == 8
        details = {"global_T_8_contrast": global_t, "global_T_redundant_MP": redundant_t,
                   "absolute_T_discrepancy": abs(global_t-redundant_t), "pair_T": pair_t,
                   "full_pair_covariance_rank": 8, "full_pair_covariance": full_cov.tolist(),
                   "independent_contrast_covariance": covariance.tolist(),
                   "shared_block_extraction_error": float(np.max(np.abs(shared-h[:, :2, :2])))}
    return global_t, np.asarray(pair_t), details


def fit_common_quotient(z, projected_y, gram, design):
    normal = np.zeros((19, 19))
    rhs = np.zeros(19)
    for c, f in enumerate(design["block_maps"]):
        normal += f.T @ np.kron(gram[c], np.eye(3)) @ f
        rhs += f.T @ (projected_y[c].T @ z[c]).reshape(9, order="F")
    theta = np.linalg.solve(normal, rhs)
    return np.stack([(f @ theta).reshape(3, 3, order="F") for f in design["block_maps"]])


def generate_observations(generator, blocks, design):
    sigma = CONFIG["output_noise_standard_deviation"]
    x_cal = generator.standard_normal((3, 30, 6))
    noise_cal = sigma*generator.standard_normal((3, 30, 6))
    x_test = generator.standard_normal((3, 20, 6))
    noise_test = sigma*generator.standard_normal((3, 20, 6))
    kernels = design["local_out"] @ blocks @ design["local_in"].transpose(0, 2, 1)
    y_cal = x_cal @ kernels.transpose(0, 2, 1) + noise_cal
    y_test = x_test @ kernels.transpose(0, 2, 1) + noise_test
    return x_cal, y_cal, x_test, y_test


def verify_design(design):
    po, pi = design["P_out"], design["P_in"]
    eo, ei = design["E_out"], design["E_in"]
    lo, li = design["local_out"], design["local_in"]
    po_plus = eo @ lo.transpose(0, 2, 1)
    pi_plus = li @ ei.transpose(0, 2, 1)
    mp_out_error = float(np.max(np.abs(po_plus-np.linalg.pinv(po))))
    mp_in_error = float(np.max(np.abs(pi_plus-np.linalg.pinv(pi))))
    assert mp_out_error < 1e-12 and mp_in_error < 1e-12
    output_projectors = po_plus @ po
    input_projectors = pi @ pi_plus
    assert all(rank_info(p)["rank"] == 3 for p in output_projectors)
    assert all(rank_info(p)["rank"] == 3 for p in input_projectors)
    common_out = design["common_out"] @ design["common_out"].T
    common_in = design["common_in"] @ design["common_in"].T
    assert rank_info(common_out)["rank"] == rank_info(common_in)["rank"] == 2
    absorption = max(float(np.max(np.abs(common_out @ output_projectors-common_out))),
                     float(np.max(np.abs(input_projectors @ common_in-common_in))))
    assert absorption < 1e-12
    portfolio = np.vstack([np.kron(p_i.T, p_o) for p_o, p_i in zip(po, pi)])
    portfolio_info = rank_info(portfolio)
    assert portfolio_info["rank"] == 19
    hidden_entries = sorted(set((i, j) for i in range(6) for j in range(6))-set(design["visible_entries"]))
    hidden_basis = np.column_stack([np.outer(design["Q_out"][:, i], design["Q_in"][:, j]).reshape(36, order="F")
                                   for i, j in hidden_entries])
    assert len(hidden_entries) == 17 and rank_info(hidden_basis)["rank"] == 17
    hidden_error = float(np.max(np.abs(portfolio @ hidden_basis)))
    assert hidden_error < 1e-12
    a, b = design["factor_A"], design["factor_B"]
    tangent_columns = []
    for i in range(6):
        for k in range(2):
            tangent_columns.append(np.outer(np.eye(6)[:, i], b[:, k]).reshape(36, order="F"))
    for j in range(6):
        for k in range(2):
            tangent_columns.append(np.outer(a[:, k], np.eye(6)[:, j]).reshape(36, order="F"))
    tangent = np.column_stack(tangent_columns)
    tangent_info, visible_tangent_info = rank_info(tangent), rank_info(portfolio @ tangent)
    assert tangent_info["rank"] == 20 and visible_tangent_info["rank"] == 16
    altered_a = a.copy()
    altered_a[:, 0] += 0.37*design["Q_out"][:, 5]
    invisible_witness = altered_a @ b.T
    witness_difference = invisible_witness-design["K_true"]
    witness_observation_error = float(np.max(np.abs(po @ witness_difference @ pi)))
    assert rank_info(invisible_witness)["rank"] == 2 and witness_observation_error < 1e-12

    check_rng = rng(CONFIG["seeds"]["deterministic_verification_probes"])
    x_cal, y_cal, x_test, y_test = generate_observations(check_rng, design["H0"], design)
    h_cal, _, z_cal, py_cal, g_cal = fit_readable_blocks(x_cal, y_cal, design)
    h_test, gi_test, z_test, _, _ = fit_readable_blocks(x_test, y_test, design)
    _, _, covariance_check = shared_statistics(h_test, gi_test, design, full_check=True)
    full_noise_covariance_error = 0.0
    for c in range(3):
        r = z_test[c] @ gi_test[c, :, :2]
        noise_operator = np.kron(r.T, lo[c, :, :2].T)
        expected = CONFIG["output_noise_standard_deviation"]**2*np.kron(gi_test[c, :2, :2], np.eye(2))
        direct = CONFIG["output_noise_standard_deviation"]**2*noise_operator @ noise_operator.T
        full_noise_covariance_error = max(full_noise_covariance_error, float(np.max(np.abs(direct-expected))))
    assert full_noise_covariance_error < 1e-12
    common_fit = fit_common_quotient(z_cal, py_cal, g_cal, design)
    full_design = np.vstack([np.kron(z_cal[c], np.eye(3)) @ design["block_maps"][c] for c in range(3)])
    targets = np.concatenate([py_cal[c].reshape(-1) for c in range(3)])
    reference_theta = np.linalg.lstsq(full_design, targets, rcond=None)[0]
    reference_fit = np.stack([(f @ reference_theta).reshape(3, 3, order="F") for f in design["block_maps"]])
    common_fit_error = float(np.linalg.norm(common_fit-reference_fit))
    assert common_fit_error < 1e-10
    kernels0 = lo @ design["H0"] @ li.transpose(0, 2, 1)
    noiseless_y = x_test @ kernels0.transpose(0, 2, 1)
    h_noiseless, _, _, _, _ = fit_readable_blocks(x_test, noiseless_y, design)
    noiseless_error = float(np.linalg.norm(h_noiseless-design["H0"]))
    assert noiseless_error < 1e-10
    lifted = eo @ design["H0"] @ ei.transpose(0, 2, 1)
    lifted_reference = output_projectors @ design["K_true"] @ input_projectors
    reconstruction_error = float(np.max(np.abs(lifted-lifted_reference)))
    mask = np.zeros((6, 6), dtype=int)
    for sector in CONFIG["readable_coordinate_sets_zero_based"]:
        mask[np.ix_(sector, sector)] += 1
    return {"P_out_pseudoinverse_error": mp_out_error, "P_in_pseudoinverse_error": mp_in_error,
            "channel_projector_rank": 3, "shared_sector_rank": 2, "projector_absorption_error": absorption,
            "partial_reconstruction_identity_error": reconstruction_error,
            "portfolio_operator": portfolio_info, "portfolio_kernel_dimension": 17,
            "hidden_coordinate_entries_zero_based": hidden_entries, "kernel_basis_observation_error": hidden_error,
            "factor_jacobian": tangent_info, "visible_factor_jacobian": visible_tangent_info,
            "generic_rank_2_nominal_dimension": 20, "generic_visible_rank_2_dimension": 16,
            "rank_2_invisible_tangent_dimension": 4,
            "rank_2_invisible_witness": {"K_altered": invisible_witness.tolist(),
                "alteration": "add 0.37 times declared sixth output basis vector to first left factor column",
                "rank": 2, "frobenius_change": float(np.linalg.norm(witness_difference)),
                "maximum_channel_kernel_change": witness_observation_error},
            "shared_covariance_check": covariance_check,
            "direct_noise_operator_covariance_error": full_noise_covariance_error,
            "common_quotient_normal_equations_vs_full_lstsq_error": common_fit_error,
            "noiseless_readable_block_recovery_error": noiseless_error,
            "observability_mask": mask.tolist(), "common_shared_block_singular_values":
                np.linalg.svd(design["H0"][0, :2, :2], compute_uv=False).tolist()}


def run_replicates(seed, design, local_delta=None):
    generator = rng(seed)
    n = 2500
    global_t, pair_t = np.empty(n), np.empty((n, 3))
    common_sse, independent_sse = np.empty(n), np.empty(n)
    realized_delta = np.empty(n) if local_delta is not None else None
    baseline_norm = float(np.linalg.norm(design["H0"]))
    max_rank_residual, min_second_singular = 0.0, float("inf")
    range_error = 0.0
    first_covariance_check = None
    for index in range(n):
        if local_delta is None:
            blocks = design["H0"]
        else:
            directions = generator.standard_normal((3, 2, 3))
            perturbation = design["perturbation_bases"] @ directions
            perturbation *= baseline_norm/np.linalg.norm(perturbation)
            blocks = design["H0"] + local_delta*perturbation
            realized_delta[index] = np.linalg.norm(blocks-design["H0"])/baseline_norm
            singulars = np.linalg.svd(blocks, compute_uv=False)
            max_rank_residual = max(max_rank_residual, float(np.max(singulars[:, 2])))
            min_second_singular = min(min_second_singular, float(np.min(singulars[:, 1])))
            if index == 0:
                kernels = design["local_out"] @ blocks @ design["local_in"].transpose(0, 2, 1)
                output_p = design["local_out"] @ design["local_out"].transpose(0, 2, 1)
                input_p = design["local_in"] @ design["local_in"].transpose(0, 2, 1)
                range_error = float(np.max(np.abs(kernels-output_p @ kernels @ input_p)))
        x_cal, y_cal, x_test, y_test = generate_observations(generator, blocks, design)
        h_cal, _, z_cal, py_cal, g_cal = fit_readable_blocks(x_cal, y_cal, design)
        h_test, ginv_test, _, _, _ = fit_readable_blocks(x_test, y_test, design)
        global_t[index], pair_t[index], check = shared_statistics(h_test, ginv_test, design, full_check=(index == 0))
        if index == 0:
            first_covariance_check = check
        h_common = fit_common_quotient(z_cal, py_cal, g_cal, design)
        independent_kernel = design["local_out"] @ h_cal @ design["local_in"].transpose(0, 2, 1)
        common_kernel = design["local_out"] @ h_common @ design["local_in"].transpose(0, 2, 1)
        sigma2 = CONFIG["output_noise_standard_deviation"]**2
        independent_sse[index] = np.sum((y_test-x_test @ independent_kernel.transpose(0, 2, 1))**2)/sigma2
        common_sse[index] = np.sum((y_test-x_test @ common_kernel.transpose(0, 2, 1))**2)/sigma2
    if local_delta is not None:
        assert np.max(np.abs(realized_delta-local_delta)) < 1e-12
        assert max_rank_residual < 1e-12 and min_second_singular > 1e-12
        assert range_error < 1e-12
    return {"global_T": global_t, "pair_T": pair_t, "common_predictive_SSE": common_sse,
            "independent_predictive_SSE": independent_sse, "realized_delta": realized_delta,
            "maximum_third_singular_value": max_rank_residual if local_delta is not None else None,
            "minimum_second_singular_value": min_second_singular if local_delta is not None else None,
            "first_replicate_range_containment_error": range_error if local_delta is not None else None,
            "first_replicate_covariance_check": first_covariance_check}


def wilson(successes, trials):
    z = CONFIG["wilson_z"]
    p = successes/trials
    denominator = 1+z*z/trials
    center = (p+z*z/(2*trials))/denominator
    half = z*np.sqrt(p*(1-p)/trials+z*z/(4*trials*trials))/denominator
    return [0.0 if successes == 0 else float(center-half), 1.0 if successes == trials else float(center+half)]


def summarize(run, threshold, evaluate_interval=True):
    values = run["global_T"]
    rejected = int(np.count_nonzero(values > threshold))
    difference = run["independent_predictive_SSE"]-run["common_predictive_SSE"]
    return {"replicates": len(values), "rejections": rejected, "empirical_rejection_rate": rejected/len(values),
            "wilson_95": wilson(rejected, len(values)) if evaluate_interval else None,
            "mean_global_T": float(np.mean(values)), "median_global_T": float(np.median(values)),
            "global_T_quantiles": {str(p): float(v) for p, v in zip([0, .025, .5, .95, .975, 1],
                np.quantile(values, [0, .025, .5, .95, .975, 1], method="linear"))},
            "mean_pair_T": np.mean(run["pair_T"], axis=0).tolist(),
            "median_pair_T": np.median(run["pair_T"], axis=0).tolist(),
            "predictive_diagnostic_only": {
                "common_model_parameters": 19, "independent_model_parameters": 27,
                "mean_standardized_SSE_common": float(np.mean(run["common_predictive_SSE"])),
                "mean_standardized_SSE_independent": float(np.mean(run["independent_predictive_SSE"])),
                "mean_SSE_independent_minus_common": float(np.mean(difference)),
                "median_SSE_independent_minus_common": float(np.median(difference)),
                "fraction_common_has_lower_SSE": float(np.mean(difference > 0)),
                "scope": "descriptive predictive comparison, not a calibrated decision rule or standard-physics comparator"}}


def serialize_run(run):
    return {key: value.tolist() if isinstance(value, np.ndarray) else value for key, value in run.items()}


def chi_square_8_quantile(probability):
    # Exact gamma survival for 8 degrees of freedom: exp(-x/2)*sum_{j=0}^3 (x/2)^j/j!.
    low, high = 0.0, 100.0
    for _ in range(80):
        value = (low+high)/2
        t = value/2
        cdf = 1-np.exp(-t)*(1+t+t*t/2+t*t*t/6)
        if cdf < probability:
            low = value
        else:
            high = value
    return (low+high)/2


def plot_results(verification, size_run, local_results, threshold, delta_min, target):
    with plt.rc_context({"font.size": 10, "axes.labelsize": 10, "pdf.fonttype": 42,
                         "axes.spines.top": False, "axes.spines.right": False}):
        fig = plt.figure(figsize=(8.2, 6.5), constrained_layout=True)
        grid = fig.add_gridspec(2, 2, height_ratios=[1, 1.05])
        mask_axis, null_axis, power_axis = fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1]), fig.add_subplot(grid[1, :])
        mask = np.asarray(verification["observability_mask"])
        colors = ListedColormap(["#F0F2F5", "#A9CDDC", "#538DA8", "#164C66"])
        norm = BoundaryNorm([-.5, .5, 1.5, 2.5, 3.5], colors.N)
        display = mask_axis.imshow(mask, cmap=colors, norm=norm)
        for i in range(6):
            for j in range(6):
                mask_axis.text(j, i, str(mask[i, j]), ha="center", va="center", fontsize=9,
                               color="white" if mask[i, j] == 3 else "#263746")
        mask_axis.set_xticks(range(6), range(1, 7))
        mask_axis.set_yticks(range(6), range(1, 7))
        mask_axis.set_xlabel("Input coordinate in fixed basis")
        mask_axis.set_ylabel("Output coordinate in fixed basis")
        mask_axis.set_title("(a) Partial coefficient coverage", loc="left", fontsize=11)
        bar = fig.colorbar(display, ax=mask_axis, orientation="horizontal", shrink=.72, pad=.06, ticks=[0, 1, 2, 3])
        bar.set_label("Number of observing channels", fontsize=9)
        values = size_run["global_T"]
        upper = float(max(np.max(values)+.5, threshold+2))
        null_axis.hist(values, bins=np.linspace(0, upper, 35), density=True, color="#2B7896", alpha=.8,
                       edgecolor="white", linewidth=.4, label="Independent null evaluation")
        curve_x = np.linspace(0, upper, 350)
        null_axis.plot(curve_x, curve_x**3*np.exp(-curve_x/2)/96, color="#575757", linewidth=1,
                       label=r"Conditional $\chi^2_8$ density")
        null_axis.axvline(threshold, color="#B25132", linestyle="--", linewidth=1.2, label="MC rejection threshold")
        null_axis.set_xlim(0, upper)
        null_axis.set_xlabel(r"Global shared-sector $T$")
        null_axis.set_ylabel("Probability density")
        null_axis.set_title("(b) Correlation-aware null calibration", loc="left", fontsize=11)
        null_axis.legend(frameon=True, facecolor="white", edgecolor="none", framealpha=.98,
                         fontsize=7.6, loc="upper right")
        deltas = np.array([record["delta"] for record in local_results])
        rates = np.array([record["summary"]["empirical_rejection_rate"] for record in local_results])
        intervals = np.array([record["summary"]["wilson_95"] for record in local_results])
        power_axis.errorbar(deltas, rates, yerr=np.stack([rates-intervals[:, 0], intervals[:, 1]-rates]),
                            color="#236A88", marker="o", markersize=4.3, linewidth=1.4, capsize=2.5,
                            label="Empirical rejection rate; pointwise Wilson 95% intervals")
        power_axis.axhline(.95, color="#B25132", linestyle="--", linewidth=1.1, label="95% planning target")
        power_axis.axhline(.05, color="#81868A", linestyle=":", linewidth=1.0, label="Nominal 5% size")
        if delta_min is not None:
            power_axis.axvline(delta_min, color="#B25132", linestyle=":", linewidth=1)
            power_axis.annotate(r"First tested $\delta$ at target: "+f"{delta_min:g}",
                                xy=(delta_min, .95), xytext=(.105, .59), fontsize=9,
                                arrowprops={"arrowstyle": "->", "color": "#B25132", "lw": .9})
        power_axis.set_xlim(-.004, .204)
        power_axis.set_ylim(-.025, 1.035)
        power_axis.set_xlabel(r"Relative RMS response displacement $\delta$")
        power_axis.set_ylabel("Rejection probability estimate")
        power_axis.set_title("(c) Local rank-two alternatives inside the frozen map ranges", loc="left", fontsize=11)
        power_axis.grid(axis="y", linewidth=.6, color="#DEE3E7")
        power_axis.set_axisbelow(True)
        power_axis.legend(frameon=True, facecolor="white", edgecolor="none", framealpha=.98,
                          fontsize=8, loc="lower right")
        fig.savefig(target, format="pdf", metadata={"Title": "V7 partial readout and local rank-two power",
            "Author": "Henry Arellano-Pena", "CreationDate": datetime(2026, 11, 25, tzinfo=timezone.utc),
            "ModDate": datetime(2026, 11, 25, tzinfo=timezone.utc)})
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--design-only", action="store_true")
    args = parser.parse_args()
    root = args.output_root.resolve()
    (root/"results").mkdir(parents=True, exist_ok=True)
    (root/"figures").mkdir(parents=True, exist_ok=True)
    declaration = {"configuration": CONFIG, "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    declaration_path = root/"results"/"bridge_partial_local_design.json"
    # This specification is written and closed before any simulation draw.
    declaration_path.write_text(json.dumps(declaration, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    if args.design_only:
        print(f"Saved pre-run computational specification: {declaration_path}")
        return
    design = construct_design()
    verification = verify_design(design)
    print("Verified partial maps, rank-2 invisible witness, covariance equivalence and estimator.", flush=True)
    calibration = run_replicates(CONFIG["seeds"]["threshold"], design)
    threshold = float(np.quantile(calibration["global_T"], .95, method="linear"))
    print(f"Calibrated threshold: {threshold:.12g}", flush=True)
    size_run = run_replicates(CONFIG["seeds"]["size"], design)
    size_summary = summarize(size_run, threshold)
    print(f"Independent null rejection rate: {size_summary['empirical_rejection_rate']:.4f}", flush=True)
    local_results = []
    for delta in DELTAS:
        run = run_replicates(CONFIG["seeds"]["local_by_delta"][str(delta)], design, local_delta=delta)
        summary = summarize(run, threshold)
        local_results.append({"delta": delta, "seed": CONFIG["seeds"]["local_by_delta"][str(delta)],
                              "summary": summary, "statistics_and_checks": serialize_run(run)})
        print(f"delta={delta:g}: {summary['rejections']}/2500 = {summary['empirical_rejection_rate']:.4f}", flush=True)
    qualifying = [record["delta"] for record in local_results if record["delta"] > 0
                  and record["summary"]["empirical_rejection_rate"] >= CONFIG["power_planning_target"]]
    delta_min = min(qualifying) if qualifying else None
    figure_path = root/"figures"/"bridge_partial_local_results.pdf"
    plot_results(verification, size_run, local_results, threshold, delta_min, figure_path)
    data_design = {key: value.tolist() if isinstance(value, np.ndarray) else value for key, value in design.items()}
    output = {"configuration": CONFIG, "script_sha256": declaration["script_sha256"],
              "software": {"python": platform.python_version(), "numpy": np.__version__,
                           "matplotlib": matplotlib.__version__, "platform": platform.platform()},
              "design": data_design, "verification": verification,
              "MC_threshold": threshold, "conditional_chi_square_df": 8,
              "conditional_chi_square_95_reference": chi_square_8_quantile(.95),
              "threshold_calibration_summary": summarize(calibration, threshold, evaluate_interval=False),
              "threshold_calibration_statistics": serialize_run(calibration),
              "null_size_summary": size_summary, "null_size_statistics": serialize_run(size_run),
              "local_power": local_results, "first_tested_delta_at_empirical_95_percent_power": delta_min,
              "planning_scope": "grid- and ensemble-specific point estimate, not worst-case, uniform, simultaneous, or interpolated",
              "figure": "figures/bridge_partial_local_results.pdf"}
    output_path = root/"results"/"bridge_partial_local_results.json"
    output_path.write_text(json.dumps(output, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(f"Saved results: {output_path}", flush=True)
    print(f"First tested delta at empirical 95% power: {delta_min}", flush=True)


if __name__ == "__main__":
    main()
