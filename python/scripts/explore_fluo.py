"""
Raw fluorescence spectra exploration — Python equivalent of treat_fluo.m

Usage: edit the PARAMETERS block below then run:
    python python/scripts/explore_fluo.py

A file picker opens (filtered on *fluo.mat). For each figure displayed:
press 's' to export PNG+SVG, 'c' to copy the figure to the clipboard.

RAW and NORMALIZED views show, for each selected file, the mean ± 1 standard
deviation across its acquisitions (one color per file, one figure per laser).
To inspect a single measurement acquisition-by-acquisition instead (e.g. to
spot probe movement between repeats), use explore_single_measurement.py.
"""

import os
import sys

import numpy as np

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PYTHON_DIR = os.path.dirname(_THIS_DIR)
for _p in (_THIS_DIR, _PYTHON_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from _plotting import apply_style, pick_files, label_from_path, enable_export, plot_acquisitions, MATLAB_COLORS
from functions.load_fluo import load_fluo

import matplotlib.pyplot as plt

# ── PARAMETERS ─────────────────────────────────────────────────────────────
DATA_DIR = r'D:\Lyon thèse\data'   # starting folder for the file picker
WL_MIN = 480
WL_MAX = 645

SHOW_RAW = True                 # mean +/- 1 sigma per file, one figure per laser
SHOW_NORMALIZED = True          # same, each acquisition normalized to area = 1 first
SHOW_COMBINED_LASERS = True     # both lasers overlaid: color = file, linestyle = laser
# ────────────────────────────────────────────────────────────────────────────


def main():
    apply_style()

    files = pick_files(DATA_DIR, '*fluo.mat', "Select one or more *fluo.mat files")
    if not files:
        print("No file selected.")
        return

    labels = [label_from_path(f) for f in files]
    colors = [MATLAB_COLORS[i % len(MATLAB_COLORS)] for i in range(len(labels))]

    all_375, all_405 = [], []
    lam = None
    for f in files:
        a375, a405, lam = load_fluo(f, WL_MIN, WL_MAX, return_all=True)
        all_375.append(a375)
        all_405.append(a405)

    lasers = [('405', "405 nm", all_405), ('375', "375 nm", all_375)]

    if SHOW_RAW:
        for key, laser_title, data_list in lasers:
            fig, ax = plt.subplots()
            for data, color, lbl in zip(data_list, colors, labels):
                plot_acquisitions(ax, lam, data, color, mode='mean', label=lbl)
            ax.set_xlim(WL_MIN, WL_MAX)
            ax.set_xlabel("Wavelength (nm)")
            ax.set_ylabel("F.I. (a.u.)")
            ax.set_title(f"{laser_title} excitation")
            ax.legend(fontsize=8)
            enable_export(fig, f'fluo_raw_{key}')

    if SHOW_NORMALIZED:
        for key, laser_title, data_list in lasers:
            fig, ax = plt.subplots()
            for data, color, lbl in zip(data_list, colors, labels):
                data = np.atleast_2d(data)
                norm_data = data / np.sum(data, axis=1, keepdims=True)
                plot_acquisitions(ax, lam, norm_data, color, mode='mean', label=lbl)
            ax.set_xlim(WL_MIN, WL_MAX)
            ax.set_xlabel("Wavelength (nm)")
            ax.set_ylabel("Normalized F.I.")
            ax.set_title(f"{laser_title} excitation")
            ax.legend(fontsize=8)
            enable_export(fig, f'fluo_norm_{key}')

    if SHOW_COMBINED_LASERS:
        fig, (ax_raw, ax_norm) = plt.subplots(1, 2, figsize=(14, 6))
        for i, lbl in enumerate(labels):
            color = colors[i]
            mean375 = np.mean(all_375[i], axis=0)
            mean405 = np.mean(all_405[i], axis=0)
            ax_raw.plot(lam, mean405, color=color, linestyle='-', label=f'{lbl} — 405 nm')
            ax_raw.plot(lam, mean375, color=color, linestyle='--', label=f'{lbl} — 375 nm')
            ax_norm.plot(lam, mean405 / np.sum(mean405), color=color, linestyle='-', label=f'{lbl} — 405 nm')
            ax_norm.plot(lam, mean375 / np.sum(mean375), color=color, linestyle='--', label=f'{lbl} — 375 nm')
        for ax, ylabel, title in [
            (ax_raw, "F.I. (a.u.)", "Raw"),
            (ax_norm, "Normalized F.I.", "Normalized"),
        ]:
            ax.set_xlim(WL_MIN, WL_MAX)
            ax.set_xlabel("Wavelength (nm)")
            ax.set_ylabel(ylabel)
            ax.set_title(title)
            ax.legend(fontsize=7)
        fig.suptitle("375 & 405 nm excitations", fontweight='bold')
        enable_export(fig, 'fluo_combined_lasers')

    plt.show()


if __name__ == '__main__':
    main()
