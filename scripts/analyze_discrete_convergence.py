"""Compare explicit 1D writing sums with the unchanged v0.1 Eq. (4)."""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from dmdlitho import normalized_dose_profile
from dmdlitho.scanning import (DMDScanParameters, build_address_cells,
                               discrete_normalized_dose_profile, profile_error_metrics)


def main() -> None:
    width_um = 10.0
    radius_um = 2.0
    # Observe six Gaussian radii beyond both feature edges. A uniform 6001-
    # point grid resolves w independently of each address pitch; it is not
    # the writing grid. The truncated far tails are negligible for this study.
    extent_um = 0.5 * width_um + 6.0 * radius_um
    positions = np.linspace(-extent_um, extent_um, 6001)
    analytical = normalized_dose_profile(positions, width_um, radius_um)
    rows = []
    selected_profiles = {}
    for period in (4, 8, 16, 32, 64):
        scan = DMDScanParameters(writing_period=period)
        pitch = scan.address_pitch_um
        centers, weights = build_address_cells(width_um, pitch)
        discrete = discrete_normalized_dose_profile(positions, width_um, radius_um, pitch)
        errors = profile_error_metrics(discrete, analytical)
        rows.append({"writing_period": period, "address_pitch_um": pitch,
                     "pitch_over_radius": pitch / radius_um,
                     "address_count": len(centers), "total_weight_um": float(weights.sum()),
                     "relative_l2_error": errors.relative_l2_error,
                     "max_absolute_error": errors.max_absolute_error})
        if period in (4, 16, 64):
            selected_profiles[period] = discrete

    # No post-hoc amplitude correction or convergence-order fit is made.
    # Periods all divide N=768. R and R/f change with A, but are excluded
    # from these normalized shapes. Stage speed is not used to create a
    # second spatial dimension or a temporal frame simulation.
    print("A | p_add (um) | p_add/w | cells | relative L2 | max absolute")
    for row in rows:
        print(f"{row['writing_period']:2d} | {row['address_pitch_um']:.6f} | "
              f"{row['pitch_over_radius']:.4f} | {row['address_count']:5d} | "
              f"{row['relative_l2_error']:.6e} | {row['max_absolute_error']:.6e}")
    print("Profile errors compare dose shape; they are not CD or EPE metrics.")

    output = Path(__file__).resolve().parents[1] / "results"
    output.mkdir(exist_ok=True)
    summary = {"width_um": width_um, "gaussian_radius_um": radius_um,
               "beam_pitch_um": 25.0, "array_columns": 768,
               "evaluation_grid": {"min_um": -extent_um, "max_um": extent_um,
                                   "point_count": len(positions)},
               "address_convention": "centered full cells with paired partial boundary cells",
               "physical_dose_scale_factored_out": True, "rows": rows}
    (output / "discrete_convergence_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    fig, ax = plt.subplots(figsize=(7.4, 4.6), layout="constrained")
    ax.plot(positions, analytical, color="black", linewidth=2.3, label="Analytical Eq. (4)")
    for period, profile in selected_profiles.items():
        pitch = next(row["address_pitch_um"] for row in rows if row["writing_period"] == period)
        ax.plot(positions, profile, linestyle="--", linewidth=1.6,
                label=f"Discrete A={period}, pitch={pitch:.3f} µm")
    ax.set(xlabel="Lateral position (µm)", ylabel="Normalized exposure q(x)",
           title="1D discrete accumulation approaches the analytical profile",
           xlim=(-12, 12))
    ax.legend(fontsize=9)
    ax.grid(alpha=0.2)
    fig.savefig(output / "discrete_vs_analytic.png", dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.4, 4.6), layout="constrained")
    pitch_ratios = [row["pitch_over_radius"] for row in rows]
    ax.loglog(pitch_ratios, [row["relative_l2_error"] for row in rows], "o-", label="Relative L2 error")
    ax.loglog(pitch_ratios, [row["max_absolute_error"] for row in rows], "s-", label="Maximum absolute error")
    ax.set(xlabel="Address pitch / Gaussian radius (p_add / w)",
           ylabel="Numerical profile error (dimensionless)",
           title="Smaller address pitch reduces dose-profile error")
    ax.legend(fontsize=10)
    ax.grid(which="both", alpha=0.2)
    fig.savefig(output / "discrete_convergence_error.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
