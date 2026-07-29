"""
Diffuse reflectance exploration (short/long wavelength) — Python equivalent
of treat_reflectance.m + the "white light" section of traitement_data.m.

Usage: edit the PARAMETERS block below then run:
    python python/scripts/explore_reflectance.py

A file picker opens (filtered on *diffuseReflectances.mat). For each figure
displayed: press 's' to export PNG+SVG, 'c' to copy the figure to the
clipboard.

Four stages are shown, in this order:
  1. Raw signal        — straight from the instrument, no spectralon division yet
  2. Reflectance       — raw signal / spectralon source spectrum
  3. Normalized        — stage 2, each spectrum normalized to area = 1
  4. Scaled reflectance — raw signal / (spectralon * scale factor), in cm^-2
Stages 2-4 require the spectralon + scale calibration files (auto-located
under data/); if they are missing only the raw signal is shown.
"""

import os
import sys

import numpy as np
from scipy.io import loadmat
from scipy.interpolate import interp1d

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PYTHON_DIR = os.path.dirname(_THIS_DIR)
for _p in (_THIS_DIR, _PYTHON_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from _plotting import apply_style, pick_files, label_from_path, enable_export, find_project_file
from functions.load_reflectance import load_reflectance
from functions.calibration_spectralon import calibration_spectralon

import matplotlib.pyplot as plt

# ── PARAMETERS ─────────────────────────────────────────────────────────────
DATA_DIR = r'D:\Lyon thèse\data'   # starting folder for the file picker
WL_MIN = 480
WL_MAX = 645

IS_SMALL_SPECTRALON = True       # divides the spectralon spectrum by 0.8
CALIBRATE = True                 # False -> only the raw signal is shown
NAME_SPECTRALON = 'spectralon_000_diffuseReflectances'  # default name if the .mat file does not specify one

SHOW_RAW_SIGNAL = True
SHOW_REFLECTANCE = True
SHOW_NORMALIZED = True
SHOW_SCALED = True
# ────────────────────────────────────────────────────────────────────────────


def _spectralon_denominators(path, spectralon_theo_path, scale_path, idx_min, idx_max, lam):
    """
    Returns (unscaled_short, unscaled_long, scaled_short, scaled_long)
    denominators for a given *diffuseReflectances.mat file, or None if its
    spectralon measurement cannot be found.
    """
    WLight = loadmat(path, squeeze_me=True, struct_as_record=False)
    directory = os.path.dirname(path)
    if 'name_spectralon' in WLight:
        spectralon_exp_path = os.path.join(directory, str(WLight['name_spectralon']))
    else:
        spectralon_exp_path = os.path.join(directory, NAME_SPECTRALON + '.mat')
    if not os.path.isfile(spectralon_exp_path):
        print(f"[Warning] Spectralon not found for {os.path.basename(path)} "
              f"({os.path.basename(spectralon_exp_path)}) — reflectance skipped.")
        return None

    SsourceShort, SsourceLong, _ = calibration_spectralon(
        spectralon_exp_path, spectralon_theo_path, IS_SMALL_SPECTRALON)
    unscaled_short = SsourceShort[idx_min:idx_max + 1]
    unscaled_long = SsourceLong[idx_min:idx_max + 1]

    scale_data = loadmat(scale_path, squeeze_me=True, struct_as_record=False)
    Scale = scale_data['Scale']
    scale_short = interp1d(np.asarray(Scale.Short.WL).flatten(), np.asarray(Scale.Short.Sy).flatten(),
                            bounds_error=False, fill_value='extrapolate')(lam)
    scale_long = interp1d(np.asarray(Scale.Long.WL).flatten(), np.asarray(Scale.Long.Sy).flatten(),
                           bounds_error=False, fill_value='extrapolate')(lam)

    return unscaled_short, unscaled_long, unscaled_short * scale_short, unscaled_long * scale_long


def main():
    apply_style()

    files = pick_files(DATA_DIR, '*diffuseReflectances.mat',
                        "Select one or more *diffuseReflectances.mat files")
    if not files:
        print("No file selected.")
        return

    spectralon_theo_path = find_project_file('Reflectance_values_array.txt')
    scale_path = find_project_file('scale_new_theory_intralipids.mat')
    can_calibrate = CALIBRATE and bool(spectralon_theo_path and scale_path)
    if CALIBRATE and not can_calibrate:
        print("[Warning] Calibration files not found — only the raw signal will be shown.")

    labels = [label_from_path(f) for f in files]
    raw_short, raw_long = [], []
    refl_short, refl_long = [], []
    scaled_short, scaled_long = [], []
    lam = None
    for f in files:
        WLShort, WLLong, lam, WLight, idx_min, idx_max, _ = load_reflectance(f, WL_MIN, WL_MAX)
        raw_short.append(WLShort)
        raw_long.append(WLLong)
        if not can_calibrate:
            continue
        denom = _spectralon_denominators(f, spectralon_theo_path, scale_path, idx_min, idx_max, lam)
        if denom is None:
            refl_short.append(None); refl_long.append(None)
            scaled_short.append(None); scaled_long.append(None)
        else:
            d_us, d_ul, d_ss, d_sl = denom
            refl_short.append(WLShort / d_us)
            refl_long.append(WLLong / d_ul)
            scaled_short.append(WLShort / d_ss)
            scaled_long.append(WLLong / d_sl)

    def _plot_pair(data_short, data_long, ylabel, key, normalize=False):
        fig, (ax_s, ax_l) = plt.subplots(1, 2, figsize=(14, 6))
        for s, l, lbl in zip(data_short, data_long, labels):
            if s is None or l is None:
                continue
            s_plot = np.asarray(s) / np.sum(s) if normalize else s
            l_plot = np.asarray(l) / np.sum(l) if normalize else l
            ax_s.plot(lam, s_plot, label=lbl)
            ax_l.plot(lam, l_plot, label=lbl)
        for ax, sub in [(ax_s, "Short"), (ax_l, "Long")]:
            ax.set_xlim(WL_MIN, WL_MAX)
            ax.set_xlabel("Wavelength (nm)")
            ax.set_ylabel(ylabel)
            ax.set_title(sub)
            ax.legend(fontsize=8)
        enable_export(fig, key)

    if SHOW_RAW_SIGNAL:
        _plot_pair(raw_short, raw_long, "Raw signal (counts)", 'refl_raw_signal')

    if can_calibrate and SHOW_REFLECTANCE:
        _plot_pair(refl_short, refl_long, "Reflectance", 'refl_reflectance')

    if can_calibrate and SHOW_NORMALIZED:
        _plot_pair(refl_short, refl_long, "Normalized reflectance", 'refl_normalized', normalize=True)

    if can_calibrate and SHOW_SCALED:
        _plot_pair(scaled_short, scaled_long, "Scaled reflectance (cm$^{-2}$)", 'refl_scaled')

    plt.show()


if __name__ == '__main__':
    main()
