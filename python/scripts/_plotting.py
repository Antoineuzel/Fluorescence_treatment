"""
Shared helpers for the exploration scripts in python/scripts/explore_*.py.

Purpose: give these scripts a consistent figure style (close to the usual
MATLAB look: thick lines, readable fonts) and a quick way to export/copy a
figure for writing up results — the PyQt5 application cannot cleanly copy its
embedded figures, so these scripts exist precisely to produce standalone
figures that are ready to paste elsewhere.
"""

import os
import io
import tkinter as tk
from tkinter import filedialog

import numpy as np
import matplotlib.pyplot as plt

# Palette identical to `matlab_colors` used in treat_fluo.m / treat_reflectance.m
MATLAB_COLORS = [
    '#0072BD', '#D95319', '#EDB120', '#7E2F8E',
    '#77AC30', '#4DBEEE', '#A2142F', '#D98080',
]

# Fit component colors/labels — identical to python/gui/main_window.py
COMPONENT_COLORS = {
    'FAD':      '#8B2FC9',
    'NADH':     '#D4A017',
    'FMN':      '#1A9970',
    'Lipo':     '#2E86DE',
    'PpIX_636': '#C0392B',
    'PpIX_620': '#E67E22',
}
COMPONENT_LABELS = {
    'FAD':      'FAD',
    'NADH':     'NADH',
    'FMN':      'protein-bound FMN',
    'Lipo':     'Lipopigments',
    'PpIX_636': 'PpIX 636 nm',
    'PpIX_620': 'PpIX 620 nm',
}


def apply_style():
    """Applies a consistent figure style (readable fonts, thick lines)."""
    plt.rcParams.update({
        'font.size': 13,
        'axes.titlesize': 14,
        'axes.titleweight': 'bold',
        'axes.labelsize': 13,
        'axes.labelweight': 'bold',
        'axes.grid': True,
        'grid.alpha': 0.3,
        'lines.linewidth': 2.5,
        'legend.fontsize': 10,
        'legend.framealpha': 0.9,
        'figure.figsize': (9, 6),
        'savefig.dpi': 300,
        'axes.prop_cycle': plt.cycler(color=MATLAB_COLORS),
    })


def label_from_path(path, strip=('fluo', 'diffuseReflectances', 'treated')):
    """Readable legend name built from a .mat file path."""
    name = os.path.splitext(os.path.basename(path))[0]
    for token in strip:
        name = name.replace(f'_{token}', '').replace(token, '')
    return name.replace('_', ' ').strip()


def find_project_file(filename, start_dir=None, max_levels=6):
    """
    Looks for `filename` by walking up from `start_dir` (defaults to this
    file's folder, i.e. python/scripts/) and scanning subfolders at each
    level — avoids hardcoding the path to data/.
    """
    if start_dir is None:
        start_dir = os.path.dirname(os.path.abspath(__file__))
    current = start_dir
    for _ in range(max_levels):
        for root, dirs, files in os.walk(current):
            dirs[:] = [d for d in dirs if d not in ('.git', '__pycache__')]
            if filename in files:
                return os.path.join(root, filename)
        parent = os.path.dirname(current)
        if parent == current:
            break
        current = parent
    return None


def pick_files(data_dir, pattern, title="Select one or more files"):
    """
    Opens a native file picker (tkinter) filtered by `pattern`
    (e.g. '*fluo.mat'), without depending on PyQt / QApplication.

    Returns
    -------
    list[str] : full paths of the selected files (can be empty)
    """
    root = tk.Tk()
    root.withdraw()
    root.attributes('-topmost', True)
    paths = filedialog.askopenfilenames(
        title=title,
        initialdir=data_dir if os.path.isdir(data_dir) else os.getcwd(),
        filetypes=[("Matching files", pattern), ("All files", "*.*")],
    )
    root.destroy()
    return list(paths)


def pick_file(data_dir, pattern, title="Select a file"):
    """
    Same as pick_files but for a single file — used by scripts that study
    one measurement at a time (e.g. explore_single_measurement.py).

    Returns
    -------
    str or None : full path of the selected file, or None if cancelled
    """
    root = tk.Tk()
    root.withdraw()
    root.attributes('-topmost', True)
    path = filedialog.askopenfilename(
        title=title,
        initialdir=data_dir if os.path.isdir(data_dir) else os.getcwd(),
        filetypes=[("Matching files", pattern), ("All files", "*.*")],
    )
    root.destroy()
    return path or None


def plot_acquisitions(ax, lam, data, color, mode='mean', label=None):
    """
    Plots either the mean ± 1 standard deviation across acquisitions
    (mode='mean') or every individual acquisition (mode='all') on `ax`.
    Mirrors the "Raw data" tab of the PyQt5 app. `data` can be 1D (a single
    spectrum) or 2D (acquisitions x wavelength).

    `label` overrides the default legend text — pass e.g. a filename when
    overlaying several measurements on the same axes, so the mean line reads
    as "<label>" and the σ band stays unlabeled (avoids one "± 1 σ" legend
    entry per file). Leave it None for the single-measurement case, where
    the default "Mean" / "± 1 σ" / "Acquisition i" labels are more useful.
    """
    data = np.atleast_2d(data)
    n = data.shape[0]
    if mode == 'mean':
        mu = np.mean(data, axis=0)
        ax.plot(lam, mu, color=color, lw=2.5, label=label or 'Mean')
        if n > 1:
            sigma = np.std(data, axis=0)
            sigma_label = '± 1 σ' if label is None else None
            ax.fill_between(lam, mu - sigma, mu + sigma, alpha=0.2, color=color, label=sigma_label)
    else:
        cmap = plt.cm.viridis(np.linspace(0.1, 0.9, max(n, 1)))
        for i in range(n):
            line_label = f'Acquisition {i + 1}' if label is None else f'{label} — acq. {i + 1}'
            ax.plot(lam, data[i], color=cmap[i], lw=1.2, alpha=0.8, label=line_label)


def _copy_figure_to_clipboard(fig):
    try:
        from PIL import Image
        import win32clipboard
    except ImportError as e:
        print(f"[Clipboard unavailable] {e} — use 's' to export to a file instead.")
        return
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=250, bbox_inches='tight')
    buf.seek(0)
    image = Image.open(buf).convert('RGB')
    out = io.BytesIO()
    image.save(out, 'BMP')
    data = out.getvalue()[14:]  # strip the BITMAPFILEHEADER -> DIB expected by the clipboard
    out.close()
    win32clipboard.OpenClipboard()
    win32clipboard.EmptyClipboard()
    win32clipboard.SetClipboardData(win32clipboard.CF_DIB, data)
    win32clipboard.CloseClipboard()
    print("Figure copied to the clipboard (Ctrl+V into Word/PowerPoint...).")


def _save_figure(fig, basename, export_dir):
    os.makedirs(export_dir, exist_ok=True)
    png_path = os.path.join(export_dir, basename + '.png')
    svg_path = os.path.join(export_dir, basename + '.svg')
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    fig.savefig(svg_path, bbox_inches='tight')
    print(f"Figure exported: {png_path}\n                 {svg_path}")


def enable_export(fig, basename, export_dir=None):
    """
    Binds two keyboard shortcuts on the matplotlib figure `fig` (call before
    `plt.show()`, once per figure):
      's' -> exports as PNG (300 dpi) + SVG into `export_dir` (defaults to
             an exports/ subfolder next to the scripts)
      'c' -> copies the rendered figure to the Windows clipboard
    """
    if export_dir is None:
        export_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'exports')

    def _on_key(event):
        if event.key == 's':
            _save_figure(fig, basename, export_dir)
        elif event.key == 'c':
            _copy_figure_to_clipboard(fig)

    fig.canvas.mpl_connect('key_press_event', _on_key)
