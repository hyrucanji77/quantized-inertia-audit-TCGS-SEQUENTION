#!/usr/bin/env python3
"""Generate the capacitor-scaling figure from the QI audit's stated model.

Model (SI units): F = 1.4e-4 * I * A / d**2.
The figure is a model calculation, not experimental data.
Output: figures/qi_capacitor_scaling.pdf, relative to the project root.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    coefficient = 1.4e-4  # N A^-1
    area_m2 = 0.01
    currents_microamp = (10.0, 30.0, 100.0)
    separations_micro_m = np.geomspace(10.0, 1000.0, 501)
    separations_m = separations_micro_m * 1e-6
    colors = ("#1f77b4", "#d95f02", "#2c7c66")

    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 10.5,
        "mathtext.fontset": "cm",
        "axes.labelsize": 11,
        "axes.titlesize": 11,
        "pdf.fonttype": 42,
    })
    fig, ax = plt.subplots(figsize=(6.8, 4.35))
    for current_microamp, color in zip(currents_microamp, colors):
        current_amp = current_microamp * 1e-6
        force_n = coefficient * current_amp * area_m2 / separations_m**2
        ax.loglog(
            separations_micro_m,
            force_n,
            linewidth=2.0,
            color=color,
            label=rf"$I={current_microamp:g}\,\mu\mathrm{{A}}$",
        )
        print(
            f"I={current_microamp:g} microamp: "
            f"F(10 micrometres)={force_n[0]:.8g} N; "
            f"F(1000 micrometres)={force_n[-1]:.8g} N"
        )

    ax.set_xlim(10, 1000)
    ax.set_ylim(1e-5, 3.0)
    ax.set_xlabel(r"Plate separation $d$ ($\mu\mathrm{m}$)")
    ax.set_ylabel(r"Predicted force $F$ (N)")
    ax.set_title(
        r"QI capacitor scaling: $F=1.4\times10^{-4}IA/d^2$, "
        r"$A=10^{-2}\,\mathrm{m}^2$"
    )
    ax.grid(True, which="major", alpha=0.23, linewidth=0.7)
    ax.grid(True, which="minor", alpha=0.09, linewidth=0.45)
    ax.legend(loc="upper right", frameon=True, framealpha=0.97)
    fig.tight_layout()

    output = Path(__file__).resolve().parents[1] / "figures" / "qi_capacitor_scaling.pdf"
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight", metadata={
        "Title": "Predicted QI capacitor force at fixed active area",
        "Author": "Henry Arellano-Pena",
        "Subject": "Stated QI model calculation; no experimental data",
        "Creator": "qi_capacitor_scaling.py",
    })
    plt.close(fig)
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
