"""
Python translation of fluo_exp_385.m

Fluorescence model: linear combination of reference spectra (FAD, NADH)
and Gaussians (FMN, lipopigments, PpIX 636 nm, PpIX 620 nm).

The Gaussian fluorophores can be replaced by a measured reference spectrum
by passing the optional `ref_spectra` parameter (dict name -> array).

Model parameters (14 parameters, Python 0-based indices):
  param[0]  : FAD amplitude
  param[1]  : NADH amplitude
  param[2]  : FMN amplitude (Gaussian or spectrum)
  param[3]  : Lipopigments amplitude (Gaussian or spectrum)
  param[4]  : PpIX 636 nm amplitude (Gaussian or spectrum)
  param[5]  : PpIX 620 nm amplitude (Gaussian or spectrum)
  param[6]  : FMN mean (nm)          — ignored if ref_spectra['FMN'] provided
  param[7]  : FMN std dev (nm)       — ignored if ref_spectra['FMN'] provided
  param[8]  : Lipopigments mean (nm) — ignored if ref_spectra['Lipo'] provided
  param[9]  : Lipopigments std dev   — ignored if ref_spectra['Lipo'] provided
  param[10] : PpIX 636 mean (nm)     — ignored if ref_spectra['PpIX_636'] provided
  param[11] : PpIX 636 std dev       — ignored if ref_spectra['PpIX_636'] provided
  param[12] : PpIX 620 mean (nm)     — ignored if ref_spectra['PpIX_620'] provided
  param[13] : PpIX 620 std dev       — ignored if ref_spectra['PpIX_620'] provided
"""

import numpy as np
from scipy.stats import norm


def fluo_exp_385(param, fluorophores, lambda_wl, ref_spectra=None):
    """
    Computes the modeled fluorescence spectrum.

    Parameters
    ----------
    param : array-like, length 14
        Model parameters (see module header)
    fluorophores : ndarray, shape (2, N)
        Normalized reference spectra (row 0: FAD, row 1: NADH)
    lambda_wl : ndarray, shape (N,)
        Wavelength vector (nm)
    ref_spectra : dict, optional
        Measured reference spectra to replace the Gaussians.
        Possible keys: 'FMN', 'Lipo', 'PpIX_636', 'PpIX_620'.
        Values: ndarray shape (N,), on the same grid as lambda_wl.

    Returns
    -------
    fluo : ndarray, shape (N,)
        Computed fluorescence spectrum
    """
    if ref_spectra is None:
        ref_spectra = {}
    param = np.asarray(param, dtype=float)
    lam   = np.asarray(lambda_wl, dtype=float)

    FMN      = (ref_spectra['FMN']
                if 'FMN' in ref_spectra
                else norm.pdf(lam, loc=param[6],  scale=abs(param[7]))  * 1000.0)
    Lipo     = (ref_spectra['Lipo']
                if 'Lipo' in ref_spectra
                else norm.pdf(lam, loc=param[8],  scale=abs(param[9]))  * 1000.0)
    PpIX_636 = (ref_spectra['PpIX_636']
                if 'PpIX_636' in ref_spectra
                else norm.pdf(lam, loc=param[10], scale=abs(param[11])) * 1000.0)
    PpIX_620 = (ref_spectra['PpIX_620']
                if 'PpIX_620' in ref_spectra
                else norm.pdf(lam, loc=param[12], scale=abs(param[13])) * 1000.0)

    fluo = (param[0] * fluorophores[0, :]
            + param[1] * fluorophores[1, :]
            + param[2] * FMN
            + param[3] * Lipo
            + param[4] * PpIX_636
            + param[5] * PpIX_620)

    return fluo


def extract_components(res_params, fluorophores, lambda_wl, ref_spectra=None):
    """
    Decomposes the spectrum into the individual contributions of each fluorophore.

    Parameters
    ----------
    res_params : ndarray, shape (14,)
        Parameters resulting from the fit (without alpha)
    fluorophores : ndarray, shape (2, N)
        FAD and NADH reference spectra
    lambda_wl : ndarray, shape (N,)
        Wavelength vector
    ref_spectra : dict, optional
        Same keys as in fluo_exp_385 — passed through as-is.

    Returns
    -------
    dict with keys: 'FAD', 'NADH', 'FMN', 'Lipo', 'PpIX_636', 'PpIX_620'
    """
    if ref_spectra is None:
        ref_spectra = {}
    p = res_params

    p_FAD     = [p[0], 0,    0,    0,    0,    0,    1,    1,    1,    1,    1,     1,     1,     1]
    p_NADH    = [0,    p[1], 0,    0,    0,    0,    1,    1,    1,    1,    1,     1,     1,     1]
    p_FMN     = [0,    0,    p[2], 0,    0,    0,    p[6], p[7], 1,    1,    1,     1,     1,     1]
    p_Lipo    = [0,    0,    0,    p[3], 0,    0,    1,    1,    p[8], p[9], 1,     1,     1,     1]
    p_PpIX636 = [0,    0,    0,    0,    p[4], 0,    1,    1,    1,    1,    p[10], p[11], 1,     1]
    p_PpIX620 = [0,    0,    0,    0,    0,    p[5], 1,    1,    1,    1,    1,     1,     p[12], p[13]]

    return {
        'FAD':      fluo_exp_385(p_FAD,     fluorophores, lambda_wl, ref_spectra),
        'NADH':     fluo_exp_385(p_NADH,    fluorophores, lambda_wl, ref_spectra),
        'FMN':      fluo_exp_385(p_FMN,     fluorophores, lambda_wl, ref_spectra),
        'Lipo':     fluo_exp_385(p_Lipo,    fluorophores, lambda_wl, ref_spectra),
        'PpIX_636': fluo_exp_385(p_PpIX636, fluorophores, lambda_wl, ref_spectra),
        'PpIX_620': fluo_exp_385(p_PpIX620, fluorophores, lambda_wl, ref_spectra),
    }


# Maps the short fluorophore names to the field-name fragments used in the
# `fluorophore` dict returned by functions.corrected_fluo.corrected_fluo_ls_wl_1
# (e.g. result['fluorophore']['flavine_385_exp']) -- same mapping as used
# inline in python/gui/main_window.py::_update_fit_tab.
FLUOROPHORE_NAME_MAP = {
    'FAD': 'flavine', 'NADH': 'NADH', 'FMN': 'gaussian',
    'Lipo': 'lipo', 'PpIX_636': 'PpIX_636', 'PpIX_620': 'PpIX_620',
}


def compute_fractions(fluorophore, lam, suffix):
    """
    Computes each active fluorophore's share of the total fluorescence, as a
    percentage, from a fit result's `fluorophore` dict (as returned by
    corrected_fluo_ls_wl_1). Mirrors the histogram computation in
    python/gui/main_window.py::_update_fit_tab.

    Parameters
    ----------
    fluorophore : dict
        `result['fluorophore']` from corrected_fluo_ls_wl_1.
    lam : ndarray
        Wavelength vector matching the component arrays.
    suffix : str
        '_385_exp' or '_405_exp', selects which excitation laser to read.

    Returns
    -------
    dict {name: percent} for every fluorophore whose component is present and
    non-zero (disabled fluorophores are simply absent from the result).
    """
    integrals = {}
    for name, mat_name in FLUOROPHORE_NAME_MAP.items():
        comp = fluorophore.get(mat_name + suffix)
        if comp is not None and np.any(comp != 0):
            integrals[name] = max(float(np.trapz(comp, lam)), 0.0)

    total = sum(integrals.values())
    if total <= 0:
        return {name: 0.0 for name in integrals}
    return {name: value / total * 100.0 for name, value in integrals.items()}
