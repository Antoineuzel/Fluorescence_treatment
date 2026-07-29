"""
Multi-point biomarker evolution — compare the fluorophore fit results across
several measurements (e.g. different timepoints, samples, or conditions).

Usage: edit the PARAMETERS block below then run:
    python python/scripts/explore_biomarkers_evolution.py

A file picker opens (filtered on *fluo.mat, multi-select). The fit
(functions.corrected_fluo.corrected_fluo_ls_wl_1 — the same one used by
explore_fit.py and the GUI's "Fit" tab) runs once per selected file, and the
results are compared across the selected points (in selection order, labelled
by filename) in three figures:

  1. NADH & FAD evolution — the two most important biomarkers, as % of total
     fluorescence, one color per fluorophore (fixed, matches the rest of the
     pipeline), solid line = 405 nm, dashed = 375 nm.
  2. Redox ratio evolution — its own figure/axis (a ratio, not a % of signal,
     so it never shares an axis with the fraction plots).
  3. All 6 fluorophores — 100%-stacked bar chart per point, one subplot per
     laser (405 nm | 375 nm), so composition and its evolution across points
     are both visible at a glance.

For each figure: press 's' to export PNG+SVG, 'c' to copy it to the clipboard.
"""

import os
import sys

import numpy as np

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PYTHON_DIR = os.path.dirname(_THIS_DIR)
for _p in (_THIS_DIR, _PYTHON_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from _plotting import (apply_style, pick_files, label_from_path,
                        enable_export, find_project_file, COMPONENT_COLORS, COMPONENT_LABELS)
from functions.load_fluo import load_fluo
from functions.corrected_fluo import corrected_fluo_ls_wl_1
from functions.fluo_model import compute_fractions, FLUOROPHORE_NAME_MAP

import matplotlib.pyplot as plt

# ── PARAMETERS ─────────────────────────────────────────────────────────────
DATA_DIR = r'D:\Lyon thèse\data'   # starting folder for the file picker
WL_MIN = 480
WL_MAX = 645

APPLY_OPTICAL_CORRECTION = False   # True -> correct by the Kim model (needs *diffuseReflectances.mat)
IS_SMALL_SPECTRALON = True
NAME_SPECTRALON = 'spectralon_000_diffuseReflectances'

ACTIVE_FLUOROPHORES = {
    'FAD': True, 'NADH': True, 'FMN': True,
    'Lipo': True, 'PpIX_636': True, 'PpIX_620': True,
}
# ────────────────────────────────────────────────────────────────────────────

color_375 = '#0072BD'
color_405 = '#D95319'


def _set_point_ticks(ax, labels):
    """X-tick labels for the point axis: full label text (so the exact
    measurement is always identifiable), small font + steep rotation so long
    names still fit without overlapping."""
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)


def _figsize_for(labels, panels=1, base_width=9, height=6):
    """
    Figure size that widens with the number of points, so tick labels have
    room. `base_width` is the minimum width of a SINGLE panel; `panels` is
    how many side-by-side subplots share this figure (the returned width
    already accounts for all of them — callers should not multiply it further).
    """
    per_panel = max(base_width, 0.9 * len(labels))
    return (per_panel * panels, height)


def _plot_nadh_fad(labels, frac_405, frac_375):
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=_figsize_for(labels))
    for name in ('NADH', 'FAD'):
        color = COMPONENT_COLORS[name]
        y405 = [f.get(name, 0.0) for f in frac_405]
        y375 = [f.get(name, 0.0) for f in frac_375]
        ax.plot(x, y405, color=color, linestyle='-', marker='o', label=f'{name} — 405 nm')
        ax.plot(x, y375, color=color, linestyle='--', marker='o', label=f'{name} — 375 nm')
    _set_point_ticks(ax, labels)
    ax.set_ylabel("Fraction of total F.I. (%)")
    ax.set_title("NADH & FAD evolution")
    ax.legend(fontsize=8)
    fig.tight_layout()
    enable_export(fig, 'biomarkers_nadh_fad')


def _plot_redox(labels, redox_405, redox_375):
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=_figsize_for(labels))
    ax.plot(x, redox_405, color=color_405, marker='o', label='405 nm')
    ax.plot(x, redox_375, color=color_375, marker='o', label='375 nm')
    _set_point_ticks(ax, labels)
    ax.set_ylabel("Redox ratio  (FAD / (FAD + NADH))")
    ax.set_title("Redox ratio evolution")
    ax.legend(fontsize=8)
    fig.tight_layout()
    enable_export(fig, 'biomarkers_redox')


def _plot_stacked(ax, labels, fractions_list, title):
    x = np.arange(len(labels))
    bottom = np.zeros(len(labels))
    for name in FLUOROPHORE_NAME_MAP:
        values = np.array([f.get(name, 0.0) for f in fractions_list])
        if not np.any(values):
            continue
        ax.bar(x, values, bottom=bottom, width=0.6,
               color=COMPONENT_COLORS[name], label=COMPONENT_LABELS[name])
        bottom += values
    _set_point_ticks(ax, labels)
    ax.set_ylabel("Fraction of total F.I. (%)")
    ax.set_ylim(0, 105)
    ax.set_title(title)
    ax.legend(fontsize=7)


def main():
    apply_style()

    fad_path = find_project_file('flavine_fluo.mat')
    nadh_path = find_project_file('NADH_fluo_4.mat')
    if not fad_path or not nadh_path:
        print("FAD/NADH reference files not found under data/fluorophores/.")
        return

    spectralon_theo_path = find_project_file('Reflectance_values_array.txt')
    scale_path = find_project_file('scale_new_theory_intralipids.mat')
    apply_correction = APPLY_OPTICAL_CORRECTION and bool(spectralon_theo_path and scale_path)
    if APPLY_OPTICAL_CORRECTION and not apply_correction:
        print("[Notice] Calibration files not found — running without optical correction.")

    fad_385, fad_405, _ = load_fluo(fad_path, WL_MIN, WL_MAX)
    nadh_385, nadh_405, _ = load_fluo(nadh_path, WL_MIN, WL_MAX)
    fluorophores_385 = np.vstack([fad_385, nadh_385])
    fluorophores_405 = np.vstack([fad_405, nadh_405])

    files = pick_files(DATA_DIR, '*fluo.mat', "Select the *fluo.mat files (data points) to compare")
    if not files:
        print("No file selected.")
        return

    labels, frac_405, frac_375, redox_405, redox_375 = [], [], [], [], []
    for f in files:
        label = label_from_path(f)
        print(f"Fitting {label}...")
        try:
            result = corrected_fluo_ls_wl_1(
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
        except Exception as e:
            print(f"[Error] {label}: {e} — skipped.")
            continue

        lam = result['lambda_fluo']
        fluo = result['fluorophore']
        labels.append(label)
        frac_405.append(compute_fractions(fluo, lam, '_405_exp'))
        frac_375.append(compute_fractions(fluo, lam, '_385_exp'))
        redox_405.append(fluo.get('redox_405', 0.0))
        redox_375.append(fluo.get('redox_385', 0.0))

    if not labels:
        print("No file could be fitted successfully.")
        return

    _plot_nadh_fad(labels, frac_405, frac_375)
    _plot_redox(labels, redox_405, redox_375)

    fig, (ax_375, ax_405) = plt.subplots(1, 2, figsize=_figsize_for(labels, panels=2, base_width=7.5))
    _plot_stacked(ax_375, labels, frac_375, "375 nm")
    _plot_stacked(ax_405, labels, frac_405, "405 nm")
    fig.suptitle("Fluorophore composition per point", fontweight='bold')
    fig.tight_layout()
    enable_export(fig, 'biomarkers_composition')

    plt.show()


if __name__ == '__main__':
    main()
