"""Reproduce Fig. 5 physical trends from Eqs. (6) and (7)."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from dmdlitho import (LithographyParameters, reconstructed_thresholds,
                      normalized_peak_dose, normalized_removed_depth,
                      clearing_write_width)


def main() -> None:
    parameters = LithographyParameters()
    threshold, clearing = reconstructed_thresholds(parameters.gamma)
    # Width zero is a plotting endpoint (no written feature), not a valid
    # feature input. Positive widths are evaluated through the physical API.
    widths = np.linspace(0, 15, 601)
    peak = np.zeros_like(widths)
    peak[1:] = normalized_peak_dose(widths[1:], parameters.gaussian_radius_um)
    removed = normalized_removed_depth(peak, threshold, clearing, parameters.gamma)
    fig, ax = plt.subplots(figsize=(7.2, 4.5), layout="constrained")
    ax.plot(widths, peak, label="Peak dose $E_{peak}/E_0$", linewidth=2)
    ax.plot(widths, removed, label="Removed depth $d_{peak}/T_0$", linewidth=2)
    ax.axvline(clearing_write_width(parameters), color="0.45", linestyle=":",
               label="Unit-dose clearing onset")
    ax.set(xlabel="Written feature width (µm)", ylabel="Normalized dose / removed depth",
           xlim=(0, 15), ylim=(0, 1.06), title="Gaussian exposure and vertical resist removal")
    ax.legend(loc="lower right", fontsize=10)
    ax.grid(alpha=0.2)
    output = Path(__file__).resolve().parents[1] / "results"
    output.mkdir(exist_ok=True)
    fig.savefig(output / "fig5_reproduced.png", dpi=200)
    plt.close(fig)
    print(f"Fig. 5: unit-dose clearing starts at {clearing_write_width(parameters):.4f} um")


if __name__ == "__main__":
    main()
