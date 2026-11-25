#!/usr/bin/env python3
"""Deterministic scalar-QI algebra and illustrative Cassini rescalings.

This calculation is not an ephemeris fit. Banik & Kroupa (2020) tested a
specific gravitational-dipole interpolator. Its total displacement norm and
radial component are retained as separate conditional sensitivity proxies.
All accelerations use SI units. Only the Python standard library is required.
"""
from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
import json
from pathlib import Path

D = Decimal
PI = D('3.141592653589793238462643383279502884197169399375105820974944592307816406286')


def excess(g: Decimal, b: Decimal, aq: Decimal) -> Decimal:
    """Physical positive excess; choose the cancellation-free quadratic form."""
    if not (g > 0 and b >= 0 and aq > 0):
        raise ValueError('Require g > 0, b >= 0, and a_Q > 0.')
    linear = g + b - aq
    discriminant = (linear * linear + 4 * aq * g).sqrt()
    if linear >= 0:
        return 2 * aq * g / (discriminant + linear)
    return (discriminant - linear) / 2


def inverse_residual(g: Decimal, delta: Decimal, aq: Decimal) -> Decimal:
    if not (g > 0 and 0 < delta <= aq):
        raise ValueError('Require g > 0 and 0 < Delta <= a_Q.')
    return (g + delta) * (aq - delta) / delta


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parents[1] / 'results' / 'qi_scalar_cassini_check.json')
    args = parser.parse_args()
    with localcontext() as context:
        context.prec = 80
        c = D('299792458')
        theta = D('8.8e26')
        aq = 2 * c * c / theta
        g_saturn = D('6.5e-5')  # Representative input, not a fitted parameter.
        reference_delta = D('.4') * D('1.2e-10')
        displacement_norm = D('7e6')  # Rounded published total in-plane norm.
        displacement_radial = D('2.13e6')  # |delta r| at the half-orbit phase.
        quoted_tracking_limit = D('160')
        epsilon_norm = reference_delta * quoted_tracking_limit / displacement_norm
        epsilon_radial = reference_delta * quoted_tracking_limit / displacement_radial

        checks = []
        maximum_force_error = D(0)
        maximum_inverse_error = D(0)
        maximum_mond_error = D(0)
        for ratio in ['1e-8', '1e-4', '.1', '1', '10', '1e4', '1e8', '1e12']:
            g = aq * D(ratio)
            for residual_ratio in ['0', '.01', '.5', '1', '2', '100', '1e10']:
                b = aq * D(residual_ratio)
                delta = excess(g, b, aq)
                a = g + delta
                reconstructed_force = a * (1 - aq / (a + b))
                force_error = abs(reconstructed_force - g) / g
                inverse_error = abs(inverse_residual(g, delta, aq) - b) / max(aq, b)
                assert a + b > aq
                assert 0 < delta <= aq * (1 + D('1e-65'))
                assert force_error < D('1e-55')
                assert inverse_error < D('1e-55')
                maximum_force_error = max(maximum_force_error, force_error)
                maximum_inverse_error = max(maximum_inverse_error, inverse_error)
                checks.append((ratio, residual_ratio))
            # Independently evaluate the simple-MOND nu map at b = a_Q.
            a_scalar = g + excess(g, aq, aq)
            a_mond = g * (D('.5') + (D('.25') + aq / g).sqrt())
            mond_error = abs(a_scalar - a_mond) / a_mond
            assert mond_error < D('1e-65')
            maximum_mond_error = max(maximum_mond_error, mond_error)

        b_values = {
            'no_residual': D('0'),
            'constant_galactic_assignment': aq,
            'unit_residual': D('1'),
            'total_norm_proxy_threshold': inverse_residual(g_saturn, epsilon_norm, aq),
            'radial_proxy_threshold': inverse_residual(g_saturn, epsilon_radial, aq),
            'assumed_1e-15_threshold': inverse_residual(g_saturn, D('1e-15'), aq),
        }
        saturn_cases = {}
        for name, b in b_values.items():
            delta = excess(g_saturn, b, aq)
            saturn_cases[name] = {
                'b_m_s2': float(b),
                'delta_m_s2': float(delta),
                'delta_over_g': float(delta / g_saturn),
                'high_acceleration_approximation_m_s2': float(aq * g_saturn / (g_saturn + b)),
                'illustrative_total_norm_displacement_m': float(displacement_norm * delta / reference_delta),
                'illustrative_radial_component_m': float(displacement_radial * delta / reference_delta),
            }

        # Wien wavelength maximum: x = 5(1 - exp(-x)), b_W = 1/x.
        x = D('5')
        for _ in range(20):
            exp_minus_x = (-x).exp()
            step = (x - 5 * (1 - exp_minus_x)) / (1 - 5 * exp_minus_x)
            x -= step
            if abs(step) < D('1e-70'):
                break
        assert abs(x - 5 * (1 - (-x).exp())) < D('1e-65')
        bw = 1 / x
        rounded_coefficient = D('.2') * PI * PI
        restored_coefficient = bw * PI * PI

        result = {
            'scope': 'Deterministic scalar algebra and conditional sensitivity proxies; no ephemeris fit, tracking-data analysis, or universal QI exclusion.',
            'sources': {
                'scalar_QI': 'https://arxiv.org/abs/1207.7007',
                'simple_MOND_mu_nu_equations_42_43': 'https://doi.org/10.12942/lrr-2012-10',
                'gravitational_dipole_Cassini_comparison': 'https://doi.org/10.1093/mnras/staa1447',
            },
            'input_parameters': {
                'c_m_s': float(c), 'theta_m': float(theta), 'a_Q_m_s2': float(aq),
                'representative_saturn_g_m_s2': float(g_saturn),
                'BK2020_reference_alpha': 0.4, 'BK2020_empirical_a0_m_s2': 1.2e-10,
                'BK2020_reference_sunward_acceleration_m_s2': float(reference_delta),
                'BK2020_rounded_total_norm_displacement_m': float(displacement_norm),
                'BK2020_radial_displacement_m': float(displacement_radial),
                'BK2020_quoted_5sigma_tracking_limit_m': float(quoted_tracking_limit),
            },
            'conditional_sensitivity_proxies': {
                'total_norm_proxy_epsilon_m_s2': float(epsilon_norm),
                'radial_component_proxy_epsilon_m_s2': float(epsilon_radial),
                'qualification': 'Total trajectory norm is not a ranging observable. Neither proxy inherits a 5-sigma QI interpretation; projections, Sun/Earth/Saturn/spacecraft responses, and fitted ephemeris parameters must be modeled.',
            },
            'saturn_scalar_cases': saturn_cases,
            'wien_roundings': {
                'wavelength_peak_x': float(x), 'b_W': float(bw),
                'two_pi_b_W': float(2 * PI * bw),
                'rounded_beta_0p2_times_pi2': float(rounded_coefficient),
                'restored_b_W_times_pi2': float(restored_coefficient),
                'rounded_beta_coefficient_percent_below_2': float(100 * (1 - rounded_coefficient / 2)),
                'restored_beta_coefficient_percent_below_2': float(100 * (1 - restored_coefficient / 2)),
                'restored_coefficient_percent_above_rounded_beta': float(100 * (restored_coefficient / rounded_coefficient - 1)),
                'below_2_percentage_denominator': 2,
                'above_rounded_beta_percentage_denominator': float(rounded_coefficient),
            },
            'verification': {
                'decimal_precision': context.prec,
                'force_branch_and_inverse_cases': len(checks),
                'simple_MOND_identity_cases': 8,
                'maximum_relative_force_error': float(maximum_force_error),
                'maximum_inverse_error_normalized_by_max_aQ_b': float(maximum_inverse_error),
                'maximum_relative_simple_MOND_error': float(maximum_mond_error),
                'all_assertions_passed': True,
            },
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(f'Wrote {args.output}; all algebra checks passed.')
    print(json.dumps({key: saturn_cases[key] for key in ['unit_residual', 'total_norm_proxy_threshold', 'radial_proxy_threshold', 'assumed_1e-15_threshold']}, indent=2))


if __name__ == '__main__':
    main()
