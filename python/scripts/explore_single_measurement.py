"""
Single-measurement acquisition check — fluorescence + reflectance.

Picks ONE *fluo.mat file (its paired *diffuseReflectances.mat is derived
automatically, same convention as the rest of the pipeline: 'fluo' ->
'diffuseReflectances' in the filename) and plots every individual
acquisition, raw, in its own figure window:
  1. Reflectance  — Short,  every acquisition
  2. Reflectance  — Long,   every acquisition
  3. Fluorescence — 375 nm, every acquisition
  4. Fluorescence — 405 nm, every acquisition

Purpose: acquisitions within a single measurement are supposed to be
repeats of the same spot — if they visibly drift apart, that usually means
the probe moved (or something else changed) between acquisitions. No
normalization or calibration is applied here on purpose, only the raw
signal, so nothing hides that drift.

Usage: edit the PARAMETERS block below then run:
    python python/scripts/explore_single_measurement.py

For each figure: press 's' to export PNG+SVG, 'c' to copy it to the clipboard.
"""

import os
import sys

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PYTHON_DIR = os.path.dirname(_THIS_DIR)
for _p in (_THIS_DIR, _PYTHON_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from _plotting import apply_style, pick_file, label_from_path, enable_export, plot_acquisitions
from functions.load_fluo import load_fluo
from functions.load_reflectance import load_reflectance

import matplotlib.pyplot as plt

# ── PARAMETERS ─────────────────────────────────────────────────────────────
DATA_DIR = r'D:\Lyon thèse\data'   # starting folder for the file picker
WL_MIN = 480
WL_MAX = 645
# ────────────────────────────────────────────────────────────────────────────

color_375 = '#0072BD'
color_405 = '#D95319'
color_short = '#0072BD'
color_long = '#D95319'


def main():
    apply_style()

    fluo_path = pick_file(DATA_DIR, '*fluo.mat', "Select one *fluo.mat measurement")
    if not fluo_path:
        print("No file selected.")
        return
    label = label_from_path(fluo_path)

    reflectance_path = os.path.join(os.path.dirname(fluo_path),
                                     os.path.basename(fluo_path).replace('fluo', 'diffuseReflectances'))
    if not os.path.isfile(reflectance_path):
        print(f"[Notice] No matching reflectance file found ({os.path.basename(reflectance_path)}) "
              "— skipping the reflectance figures.")
    else:
        all_short, all_long, lam_refl, _, _, _, _ = load_reflectance(
            reflectance_path, WL_MIN, WL_MAX, return_all=True)

        for data, color, key, channel in [(all_short, color_short, 'short', "Short"),
                                           (all_long, color_long, 'long', "Long")]:
            fig, ax = plt.subplots()
            plot_acquisitions(ax, lam_refl, data, color, mode='all')
            ax.set_xlim(WL_MIN, WL_MAX)
            ax.set_xlabel("Wavelength (nm)")
            ax.set_ylabel("Raw signal (counts)")
            ax.set_title(f"{label} — {channel} — all acquisitions")
            ax.legend(fontsize=8)
            enable_export(fig, f'single_{label}_refl_{key}')

    all_375, all_405, lam_fluo = load_fluo(fluo_path, WL_MIN, WL_MAX, return_all=True)

    for data, color, key, laser in [(all_375, color_375, '375', "375 nm"),
                                     (all_405, color_405, '405', "405 nm")]:
        fig, ax = plt.subplots()
        plot_acquisitions(ax, lam_fluo, data, color, mode='all')
        ax.set_xlim(WL_MIN, WL_MAX)
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("F.I. (a.u.)")
        ax.set_title(f"{label} — {laser} — all acquisitions")
        ax.legend(fontsize=8)
        enable_export(fig, f'single_{label}_fluo_{key}')

    plt.show()


if __name__ == '__main__':
    main()
