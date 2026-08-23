#!/usr/bin/env python3
"""
Regenerates Figure 1 of the QI audit as a two-panel figure.

The original single-panel version showed only Reading E (explicit/kinematic),
in which m_I/m = 1 - a0/|a| goes negative below the threshold. Section 2.2 of
the amended manuscript shows that the scalar law is implicit in a, and that
closing it on the realized acceleration gives

    a = F/m + a0 ,      m_I/m = F/(F + m a0)  in  [0, 1)

so the negative branch is never reached. The caption now refers to both
readings, so the figure should show both.

Output: figures/qi_scalar_response.pdf
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, (axE, axR) = plt.subplots(1, 2, figsize=(10.6, 4.1))

# ---------------------------------------------------------------- Reading E
x = np.logspace(-1, 6, 4000)          # |a| / a0, externally supplied
eta = 1.0 - 1.0 / x
pos, neg = x >= 1.0, x <= 1.0

axE.semilogx(x[pos], eta[pos], lw=2.0, color="#1f77b4",
             label=r"$\eta_{\rm QI}=1-a_0/|a|$ (positive branch)")
axE.semilogx(x[neg], eta[neg], lw=2.0, ls="--", color="#ff7f0e",
             label=r"naive continuation below $a_0$")
axE.axhline(0.0, lw=1.0, color="#2b6ca3")
axE.axvline(1.0, lw=0.9, ls=":", color="0.35")
axE.set_xlim(1e-1, 1e6)
axE.set_ylim(-4, 1.15)
axE.set_xlabel(r"Acceleration ratio $|a|/a_0^{\rm QI}$  (supplied externally)")
axE.set_ylabel(r"Effective inertial factor $m_{\rm I}/m$")
axE.set_title("Reading E (explicit, kinematic)", fontsize=11)
axE.legend(loc="lower right", fontsize=8.5, framealpha=0.95)
axE.grid(True, which="both", alpha=0.16)
axE.text(0.16, -3.35, "negative branch\n(unphysical)", fontsize=8.5,
         color="#b35300", va="center")

# ---------------------------------------------------------------- Reading R
#   realized acceleration a = F/m + a0  ->  a/a0 = 1 + phi,  phi := F/(m a0)
phi = np.logspace(-3, 6, 4000)
a_over_a0 = 1.0 + phi
eta_R = phi / (1.0 + phi)             # = F / (F + m a0)

axR.semilogx(a_over_a0, eta_R, lw=2.0, color="#1f77b4",
             label=r"$m_{\rm I}/m=F/(F+m\,a_0^{\rm QI})$")
axR.axhline(0.0, lw=1.0, color="#2b6ca3")
axR.axvline(1.0, lw=0.9, ls=":", color="0.35")
axR.axvspan(1e-1, 1.0, color="0.85", alpha=0.7, zorder=0)
axR.set_xlim(1e-1, 1e6)
axR.set_ylim(-4, 1.15)
axR.set_xlabel(r"Realized ratio $|a|/a_0^{\rm QI}=1+F/(m\,a_0^{\rm QI})$")
axR.set_title("Reading R (implicit, realized)", fontsize=11)
axR.legend(loc="lower right", fontsize=8.5, framealpha=0.95)
axR.grid(True, which="both", alpha=0.16)
axR.text(0.135, -1.6, "dynamically\ninaccessible\nfor $F\\geq 0$",
         fontsize=8.5, color="0.25", va="center")
axR.annotate(r"$F\to 0:\ a\to a_0^{\rm QI}$", xy=(1.0, 0.0),
             xytext=(6.0, -1.15), fontsize=8.5,
             arrowprops=dict(arrowstyle="->", lw=0.9, color="0.3"))

fig.suptitle("Scalar QI response under the two closures of the implicit law",
             fontsize=12.5, y=0.99)
fig.tight_layout(rect=[0, 0, 1, 0.955])
fig.savefig("figures/qi_scalar_response.pdf", bbox_inches="tight")
print("wrote figures/qi_scalar_response.pdf")
