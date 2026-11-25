# Figures

All five figures for the V8 computational release are included. They are generated
from the stated scalar models or synthetic experiments and contain no new
observational dataset.

| File | Content | Produced by |
|---|---|---|
| [`qi_scalar_response.pdf`](qi_scalar_response.pdf) | Explicit kinematic and implicit realized readings of the scalar QI law | `scripts/qi_scalar_response.py` |
| [`qi_capacitor_scaling.pdf`](qi_capacitor_scaling.pdf) | Capacitor force scaling under the stated coefficient convention | `scripts/qi_capacitor_scaling.py` |
| [`bridge_synthetic_residuals.pdf`](bridge_synthetic_residuals.pdf) | Original three-channel residual test, with a version-neutral title | `code/bridge_synthetic_demo.py` |
| [`bridge_partial_local_results.pdf`](bridge_partial_local_results.pdf) | Partial coefficient coverage, correlated shared-sector statistic, and the complete local power grid | `code/bridge_partial_local_demo.py` |
| [`bridge_four_channel_results.pdf`](bridge_four_channel_results.pdf) | Four-channel residual test against independent rank-two blocks within the frozen ranges | `code/bridge_four_channel_demo.py` |

Run the scripts from the repository root using the commands in the
[main README](../README.md#reproducing-the-computations). Regeneration overwrites
the corresponding PDFs unless a supported alternative output path is supplied.

The three-channel and partial-readability designs fail the conditional capacity
ceiling `rho_max = 0.25`; their plotted residual behavior does not override that
decision. The four-channel design passes with `8/36`, but its plot measures power
against a separated alternative ensemble. The local power curve belongs to the
partial-readability design, not the four-channel control.

The earlier manuscript sources retained elsewhere in this repository have not
been updated to V8. Figure numbering in those historical sources may differ from
the current companion manuscripts; use the filenames and experiment descriptions
above to identify the computational outputs.

Author-generated figures are released under CC BY 4.0; see [`LICENSE`](../LICENSE).
