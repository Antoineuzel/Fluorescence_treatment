# Fluorescence_treatment

Analysis pipeline for tissue optics spectroscopy data (PhD thesis project):
autofluorescence spectra (375 / 405 nm excitation) and diffuse reflectance
spectra (short / long wavelength), with a fluorophore-decomposition fit
(FAD, NADH, protein-bound FMN, lipopigments, PpIX 636/620 nm) corrected by
the tissue's optical properties (Kim diffusion model).

The same pipeline exists twice: as the original **MATLAB** code, and as a
**Python** translation with an added GUI. Both stay functionally equivalent
and are kept in sync.

## Repository layout

```
matlab/                 MATLAB code
  functions/            Core, reusable functions (loading, calibration, fit model)
  treat_fluo.m          Explore raw fluorescence spectra across files
  treat_reflectance.m   Explore diffuse reflectance across files (4 stages, see below)
  explore_single_measurement.m
                        Check acquisition-to-acquisition repeatability of ONE measurement
  fit_fluo_alpha_fitted_LS_WL_1_unique_acq_function.m
                        Fluorophore decomposition fit for a single acquisition
  traitement_data.m     Legacy script, superseded by the files above (kept for reference)

python/
  gui/main_window.py    Interactive PyQt5 app (see "GUI" below)
  main.py               Entry point for the GUI
  functions/            Python translation of matlab/functions/
  scripts/              Lightweight parameter-block scripts, the Python
                        counterparts of the matlab/*.m exploration scripts
                        (see "Exploration scripts" below)

data/                   Reference data shipped with the pipeline (not experimental data)
  spectralon/           Theoretical spectralon reflectance standard
  fluorophores/         Reference FAD / NADH emission spectra
  probe_calibration/    Probe scale/calibration factors (Scale.Short / Scale.Long)
  absorbers/            Absorption spectra used by the optical-property fit
```

Actual experimental data (`*fluo.mat`, `*diffuseReflectances.mat`,
`*_diffuseReflectances.mat` spectralon measurements) is **not** part of this
repository — every script/tool below asks you to pick those files from
wherever they are stored on disk.

## Why both a GUI and standalone scripts?

- The **GUI** (`python/gui/main_window.py`) is the "production" tool: pick a
  `*fluo.mat` file, and it auto-fits and displays raw spectra, reflectance,
  and the fluorophore decomposition together, live, as you change filters
  (λ range, active fluorophores, optical correction).
- The **exploration scripts** (`matlab/*.m` and `python/scripts/explore_*.py`)
  exist because the GUI's embedded plots can't be copied or exported cleanly
  for a manuscript. They are simple, parameter-block scripts (edit the
  values at the top, run, get a plain matplotlib/MATLAB figure window) built
  for producing clean, standalone, exportable figures rather than for
  interactive filtering. Every figure they produce supports two keyboard
  shortcuts:
  - **`s`** → export the figure as PNG (300 dpi) + a vector file (SVG in
    Python, PDF in MATLAB)
  - **`c`** → copy the figure straight to the clipboard (paste into
    Word/PowerPoint/LaTeX)

## GUI (Python)

```
cd "d:\Lyon thèse\soft\Fluorescence_treatment"
python python/main.py
```

Requires: `numpy`, `scipy`, `matplotlib`, `PyQt5`.

Three tabs, all driven by one selected `*fluo.mat` file:
- **Raw data** — 375/405 nm spectra, mean ± 1σ or every acquisition.
- **Reflectance** — short/long wavelength reflectance (auto-derived from the
  paired `*diffuseReflectances.mat` file), visualization only, with an
  optional spectralon calibration toggle.
- **Fluorescence fit** — live fluorophore decomposition, with checkboxes to
  include/exclude each fluorophore, an optical-correction toggle, and the
  resulting redox ratio (FAD / (FAD + NADH)) at 375 and 405 nm.

## Exploration scripts

Each script has a `# PARAMETERS` / `%% PARAMETERS` block at the top (data
folder, λ range, which stages to show, active fluorophores...) — edit it,
then run the script. A file picker opens for selecting the input file(s).

| Purpose | MATLAB | Python |
|---|---|---|
| Raw fluorescence spectra, mean ± 1σ per file, normalized views, 375/405 combined | `matlab/treat_fluo.m` | `python/scripts/explore_fluo.py` |
| Diffuse reflectance, 4 stages (see below) | `matlab/treat_reflectance.m` | `python/scripts/explore_reflectance.py` |
| Fluorophore decomposition fit, with/without optical correction | `matlab/fit_fluo_alpha_fitted_LS_WL_1_unique_acq_function.m` | `python/scripts/explore_fit.py` |
| Acquisition-to-acquisition check for a single measurement (spot probe movement between repeats) | `matlab/explore_single_measurement.m` | `python/scripts/explore_single_measurement.py` |

Run a Python script with, e.g.:
```
python python/scripts/explore_fluo.py
```
Requires the same packages as the GUI, plus `pywin32` and `Pillow` for the
clipboard-copy shortcut (falls back to file export only if either is
missing). MATLAB scripts require R2020a+ (`copygraphics` / `exportgraphics`).

### Reflectance stages

`treat_reflectance.m` / `explore_reflectance.py` show up to four stages, in
this order, so raw counts are never confused with an actual reflectance:

1. **Raw signal (counts)** — straight from the instrument, no spectralon
   division yet.
2. **Reflectance** — raw signal ÷ spectralon source spectrum (dimensionless).
3. **Normalized reflectance** — stage 2, each spectrum normalized to area = 1.
4. **Scaled reflectance (cm⁻²)** — raw signal ÷ (spectralon × scale factor),
   shown last since it's the one carrying physical units.

Stages 2-4 need the spectralon + scale calibration files (auto-located
under `data/`); if they're missing, only the raw signal is shown.

## Column/file conventions (raw acquisition data)

- `donnees_acq_brut.data` column 4 (MATLAB) / column 3 (Python, 0-indexed):
  385 nm laser (fluorescence) or short-wavelength LED (reflectance)
- column 5 / 4: 405 nm laser (fluorescence) or long-wavelength LED (reflectance)
- column 7 / 6: optical power
- column 8 / 7: acquisition time (ms)
- A reflectance file is the fluorescence file's name with `fluo` replaced by
  `diffuseReflectances` (e.g. `sample_001_fluo.mat` →
  `sample_001_diffuseReflectances.mat`) — every tool derives it automatically.
