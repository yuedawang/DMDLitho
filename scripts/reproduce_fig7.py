"""Recompute thresholds for each contrast; report explicit sampling inference."""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from dmdlitho import (LithographyParameters, reconstructed_thresholds,
                      resolution_limits, solve_width_correction, solve_dose_correction)


def main() -> None:
    contrasts = np.array([0.46, 0.47, 0.48, 0.49, 0.50])
    sampling_pitch_um = 0.1
    rows = []
    for gamma in contrasts:
        parameters = LithographyParameters(gamma=float(gamma))
        # Every gamma has a new Eth/E0 and Ecl/E0 under the documented inference.
        threshold, clearing = reconstructed_thresholds(float(gamma))
        limits = resolution_limits(parameters)
        # The paper reports widths to 0.1 µm but does not state a sweep pitch.
        # This grid is an explicit reproduction choice. Ceiling selects the
        # first feasible design/target width; no published limits are hard-coded.
        sampled = sampling_pitch_um * np.ceil(np.array([
            limits.uncorrected_design_um, limits.width_corrected_um,
            limits.dose_corrected_um]) / sampling_pitch_um)
        width = solve_width_correction(float(sampled[1]), parameters)
        dose = solve_dose_correction(float(sampled[2]), parameters)
        rows.append({"gamma": float(gamma), "threshold_ratio": threshold,
                     "clearing_ratio": clearing, "continuous_limits_um": vars(limits),
                     "uncorrected_design_um": float(sampled[0]),
                     "width_target_um": float(sampled[1]), "dose_target_um": float(sampled[2]),
                     "width_reinforcement_percent": width.reinforcement_percent,
                     "dose_reinforcement_percent": dose.reinforcement_percent})
    output = Path(__file__).resolve().parents[1] / "results"
    output.mkdir(exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 4.5), layout="constrained")
    for key, label, marker in [("uncorrected_design_um", "Uncorrected design", "o"),
                               ("width_target_um", "Width-corrected target", "s"),
                               ("dose_target_um", "Dose-corrected target", "^")]:
        ax.plot(contrasts, [row[key] for row in rows], marker=marker, label=label)
    ax.set(xlabel="Photoresist contrast γ", ylabel="Minimum feasible design / target width (µm)",
           title="Resolution versus contrast (0.1 µm design grid)")
    ax.set_xticks(contrasts)
    ax.legend(fontsize=10)
    ax.grid(alpha=0.2)
    fig.savefig(output / "fig7a_reproduced.png", dpi=200)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(7.2, 4.5), layout="constrained")
    ax.plot(contrasts, [row["width_reinforcement_percent"] for row in rows],
            marker="s", label="Width reinforcement")
    ax.plot(contrasts, [row["dose_reinforcement_percent"] for row in rows],
            marker="^", label="Dose reinforcement")
    ax.set(xlabel="Photoresist contrast γ", ylabel="Required reinforcement (%)",
           title="Reinforcement at first feasible grid target")
    ax.set_xticks(contrasts)
    ax.legend(fontsize=10)
    ax.grid(alpha=0.2)
    fig.savefig(output / "fig7b_reproduced.png", dpi=200)
    plt.close(fig)
    (output / "fig7_summary.json").write_text(
        json.dumps({"sampling_pitch_um": sampling_pitch_um, "rows": rows}, indent=2) + "\n",
        encoding="utf-8")
    print("gamma | uncorrected / width / dose targets (um) | width / dose (%)")
    for row in rows:
        print(f"{row['gamma']:.2f} | {row['uncorrected_design_um']:.1f} / "
              f"{row['width_target_um']:.1f} / {row['dose_target_um']:.1f} | "
              f"{row['width_reinforcement_percent']:.2f} / {row['dose_reinforcement_percent']:.2f}")


if __name__ == "__main__":
    main()
