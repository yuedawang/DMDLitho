"""Compute Fig. 6 curves independently of the published regression targets."""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from dmdlitho import (LithographyParameters, compute_cd, matching_index,
                      solve_width_correction, solve_dose_correction,
                      resolution_limits)


def main() -> None:
    parameters = LithographyParameters()
    limits = resolution_limits(parameters)
    widths = np.linspace(5, 12, 351)
    matching = np.full_like(widths, np.nan)
    width_reinforcement = np.full_like(widths, np.nan)
    dose_reinforcement = np.full_like(widths, np.nan)
    for i, width in enumerate(widths):
        # Absent crossings remain gaps, rather than fabricated plotted CDs.
        if width >= limits.uncorrected_design_um:
            matching[i] = matching_index(width, compute_cd(width, parameters).cd_um)
        if width >= limits.width_corrected_um:
            width_reinforcement[i] = solve_width_correction(width, parameters).reinforcement_percent
        if width >= limits.dose_corrected_um:
            dose_reinforcement[i] = solve_dose_correction(width, parameters).reinforcement_percent
    fig, ax = plt.subplots(figsize=(7.2, 4.5), layout="constrained")
    ax.plot(widths, matching, label="Uncorrected matching index", linewidth=2)
    ax.plot(widths, width_reinforcement, label="Width reinforcement", linewidth=2)
    ax.plot(widths, dose_reinforcement, label="Dose reinforcement", linewidth=2)
    ax.axhline(100, color="0.5", linestyle=":")
    ax.set(xlabel="Design / target feature width (µm)",
           ylabel="Matching index / reinforcement (%)", xlim=(5, 12),
           title="Width and dose pre-compensation")
    ax.legend(fontsize=10)
    ax.grid(alpha=0.2)

    # These coordinates are published comparison locations, not curve data
    # or answers. Both inverse problems still use the package's physical model.
    cd = compute_cd(6.6, parameters)
    width_result = solve_width_correction(5.9, parameters)
    dose_result = solve_dose_correction(5.4, parameters)
    summary = {
        "threshold_mode": "reconstructed",
        "matching_at_6_6_percent": matching_index(6.6, cd.cd_um),
        "printed_cd_at_6_6_um": cd.cd_um,
        "write_width_for_5_9_um": width_result.write_width_um,
        "width_reinforcement_at_5_9_percent": width_result.reinforcement_percent,
        "dose_scale_for_5_4": dose_result.dose_scale,
        "dose_reinforcement_at_5_4_percent": dose_result.reinforcement_percent,
        "continuous_limits_um": vars(limits),
    }
    output = Path(__file__).resolve().parents[1] / "results"
    output.mkdir(exist_ok=True)
    fig.savefig(output / "fig6_reproduced.png", dpi=200)
    plt.close(fig)
    (output / "fig6_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print("DMDLitho Fig. 6 reproduction\n--------------------------------")
    print(f"Matching index at 6.6 um: {summary['matching_at_6_6_percent']:.2f}%")
    print(f"Width correction at 5.9 um: {width_result.reinforcement_percent:.2f}% "
          f"(write {width_result.write_width_um:.4f} um)")
    print(f"Dose correction at 5.4 um: {dose_result.reinforcement_percent:.2f}% "
          f"(scale {dose_result.dose_scale:.5f})")
    print("Continuous clearing limits:", json.dumps(vars(limits)))


if __name__ == "__main__":
    main()
