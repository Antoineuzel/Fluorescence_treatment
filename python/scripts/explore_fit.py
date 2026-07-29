"""
Fluorescence fit exploration (fluorophore decomposition) — Python equivalent
of fit_fluo_alpha_fitted_LS_WL_1_unique_acq_function.m

Usage: edit the PARAMETERS block below then run:
    python python/scripts/explore_fit.py

A file picker opens (filtered on *fluo.mat). One or several files can be
processed in a row. For each file, the fit is run WITHOUT optical correction
and, whenever a matching *diffuseReflectances.mat + calibration files are
available, WITH optical correction too — both are shown so the effect of the
correction can be compared. Each excitation wavelength (375 nm / 405 nm) is
drawn in its own separate figure window, so a single plot can be copied on
its own (e.g. for a paper) without dragging along the other one.

For each figure: press 's' to export PNG+SVG, 'c' to copy it to the
clipboard.

This script does not reimplement the fit model: it calls
functions.corrected_fluo.corrected_fluo_ls_wl_1, the same function used by
the "Fit" tab of the PyQt5 application.
"""

import os
import sys

import numpy as np

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PYTHON_DIR = os.path.dirname(_THIS_DIR)
for _p in (_THIS_DIR, _PYTHON_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from _plotting import (apply_style, pick_files, label_from_path, enable_export,
                        find_project_file, COMPONENT_COLORS, COMPONENT_LABELS)
from functions.load_fluo import load_fluo
from functions.corrected_fluo import corrected_fluo_ls_wl_1

import matplotlib.pyplot as plt

# ── PARAMETERS ─────────────────────────────────────────────────────────────
DATA_DIR = r'D:\Lyon thèse\data'   # starting folder for the file picker
WL_MIN = 480
WL_MAX = 645

IS_SMALL_SPECTRALON = True
NAME_SPECTRALON = 'spectralon_000_diffuseReflectances'

# Fluorophores included in the fit — set to False to exclude one
ACTIVE_FLUOROPHORES = {
    'FAD': True, 'NADH': True, 'FMN': True,
    'Lipo': True, 'PpIX_636': True, 'PpIX_620': True,
}
# ────────────────────────────────────────────────────────────────────────────

_LASER_INFO = {
    '375': ('_385_exp', 'S385total', 'S385_corrected', '375 nm'),
    '405': ('_405_exp', 'S405total', 'S405_corrected', '405 nm'),
}
_NAME_MAP = {'NADH': 'NADH', 'FAD': 'flavine', 'FMN': 'gaussian',
             'Lipo': 'lipo', 'PpIX_620': 'PpIX_620', 'PpIX_636': 'PpIX_636'}


def _plot_single(result, mode_label, mode_key, laser_key, file_label):
    lam = result['lambda_fluo']
    fluo = result['fluorophore']
    suffix, raw_key, corr_key, laser_title = _LASER_INFO[laser_key]

    fig, ax = plt.subplots(figsize=(9, 6))
    S_raw = result[raw_key]
    S_corr = result.get(corr_key, S_raw)
    norm = max(np.sum(S_corr), 1e-12)

    if mode_key == 'no_correction':
        ax.plot(lam, S_raw / norm, '--', color='#808B96', alpha=0.8, label='Raw data (normalized)')
    else:
        ax.plot(lam, S_corr / norm, color='#2C3E50', lw=2.5, label='Corrected signal', zorder=5)

    sum_fit = sum(fluo.get(_NAME_MAP[k] + suffix, np.zeros_like(lam)) for k in _NAME_MAP)
    ax.plot(lam, sum_fit / norm, color='#2471A3', lw=2.5, alpha=0.9, label='Sum of fit', zorder=4)

    for k, mat_name in _NAME_MAP.items():
        comp = fluo.get(mat_name + suffix)
        if comp is not None and np.any(comp != 0):
            ax.plot(lam, comp / norm, '--', color=COMPONENT_COLORS[k],
                    label=COMPONENT_LABELS[k], zorder=3)

    ax.set_xlim(WL_MIN, WL_MAX)
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Normalized intensity (a.u.)")
    r385, r405 = fluo.get('redox_385'), fluo.get('redox_405')
    ax.set_title(f"{file_label} — {laser_title} — {mode_label}\n"
                 f"redox 375 nm = {r385:.3f}   |   redox 405 nm = {r405:.3f}")
    ax.legend(fontsize=8)
    enable_export(fig, f'fit_{file_label}_{laser_key}_{mode_key}')


def _run_fit(f, fluorophores_385, fluorophores_405, apply_correction,
             spectralon_theo_path, scale_path):
    return corrected_fluo_ls_wl_1(
        path=os.path.dirname(f), file=os.path.basename(f),
        name_spectralon=NAME_SPECTRALON,
        is_small_spectralon=IS_SMALL_SPECTRALON,
        min_wl_reflectance=WL_MIN, max_wl_reflectance=WL_MAX,
        min_wl_fluo=WL_MIN, max_wl_fluo=WL_MAX,
        fluorophores_385=fluorophores_385, fluorophores_405=fluorophores_405,
        apply_optical_correction=apply_correction,
        active_fluorophores=ACTIVE_FLUOROPHORES,
        spectralon_theo_path=spectralon_theo_path, scale_path=scale_path,
    )


def main():
    apply_style()

    fad_path = find_project_file('flavine_fluo.mat')
    nadh_path = find_project_file('NADH_fluo_4.mat')
    if not fad_path or not nadh_path:
        print("FAD/NADH reference files not found under data/fluorophores/.")
        return

    spectralon_theo_path = find_project_file('Reflectance_values_array.txt')
    scale_path = find_project_file('scale_new_theory_intralipids.mat')
    calibration_available = bool(spectralon_theo_path and scale_path)
    if not calibration_available:
        print("[Notice] Calibration files not found — only the 'no correction' fit will be shown.")

    fad_385, fad_405, _ = load_fluo(fad_path, WL_MIN, WL_MAX)
    nadh_385, nadh_405, _ = load_fluo(nadh_path, WL_MIN, WL_MAX)
    fluorophores_385 = np.vstack([fad_385, nadh_385])
    fluorophores_405 = np.vstack([fad_405, nadh_405])

    files = pick_files(DATA_DIR, '*fluo.mat', "Select one or more *fluo.mat files to fit")
    if not files:
        print("No file selected.")
        return

    for f in files:
        label = label_from_path(f)
        print(f"Fitting {label}...")

        try:
            result_no_corr = _run_fit(f, fluorophores_385, fluorophores_405, False,
                                       spectralon_theo_path, scale_path)
        except Exception as e:
            print(f"[Error] {label} (no correction): {e}")
            continue
        for laser_key in ('375', '405'):
            _plot_single(result_no_corr, 'no correction', 'no_correction', laser_key, label)

        if not calibration_available:
            continue
        try:
            result_corr = _run_fit(f, fluorophores_385, fluorophores_405, True,
                                    spectralon_theo_path, scale_path)
        except FileNotFoundError as e:
            print(f"[Notice] {label}: {e} — skipping the optical-correction fit for this file.")
            continue
        except Exception as e:
            print(f"[Error] {label} (with correction): {e}")
            continue
        for laser_key in ('375', '405'):
            _plot_single(result_corr, 'with correction', 'with_correction', laser_key, label)

    plt.show()


if __name__ == '__main__':
    main()
