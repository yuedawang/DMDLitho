"""Reconstruct Section 4.2 scan parameters without simulating its corners."""

import json
from pathlib import Path

from dmdlitho.scanning import DMDScanParameters


def main() -> None:
    scan = DMDScanParameters()
    summary = {
        "input_parameters": vars(scan),
        "tilt_angle_deg": scan.tilt_angle_deg,
        "address_pitch_um": scan.address_pitch_um,
        "moving_pitch_um": scan.moving_pitch_um,
        "relative_pitch_mismatch": scan.relative_pitch_mismatch,
        "redundancy": scan.redundancy,
        "repeated_writing_time_s": scan.repeated_writing_time_s,
    }
    print("Hur & Seo scan-geometry reconstruction (Section 4.2)")
    print(f"Tilt: {scan.tilt_angle_deg:.5f} degrees")
    print(f"Address pitch: {scan.address_pitch_um:.6f} um")
    print(f"Moving pitch: {scan.moving_pitch_um:.6f} um")
    print(f"Relative synchronization mismatch: {100 * scan.relative_pitch_mismatch:.4f}%")
    print(f"Redundancy: {scan.redundancy}")
    print(f"Physical repeated-writing time R/f: {scan.repeated_writing_time_s:.7f} s")
    print("The paper reports both pitches as about 1.5 um; values above retain equation precision.")
    output = Path(__file__).resolve().parents[1] / "results"
    output.mkdir(exist_ok=True)
    (output / "scanning_geometry_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
