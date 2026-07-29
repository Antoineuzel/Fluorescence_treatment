"""
Python translation of load_reflectance.m
Loads and processes diffuse reflectance data from a *diffuseReflectances.mat file.

Same column convention as load_fluo:
  - donnees_acq_brut.data(:, 4:6) [MATLAB] -> columns 3,4,5 [Python]
  - WLallummee col 0 (column 3) = short wavelength LED (Short)
  - WLallummee col 1 (column 4) = long wavelength LED (Long)
"""

import numpy as np
from scipy.io import loadmat


def _find_wl_index(lambda_arr, wl_target):
    return int(np.argmin((lambda_arr - wl_target) ** 2))


def load_reflectance(path, wl_min=480, wl_max=645, return_all=False):
    """
    Loads reflectance data from a *diffuseReflectances.mat file.

    Translated from load_reflectance.m (Antoine Uzel).

    Parameters
    ----------
    path : str
        Full path to the *diffuseReflectances.mat file
    wl_min : float
        Minimum wavelength (nm)
    wl_max : float
        Maximum wavelength (nm)
    return_all : bool
        If True, returns all acquisitions (2D) instead of the average (1D)

    Returns
    -------
    WLShort : ndarray
        Short wavelength reflectance — (N,) if return_all=False, (M, N) otherwise
    WLLong : ndarray
        Long wavelength reflectance — (N,) if return_all=False, (M, N) otherwise
    lambda_wl : ndarray, shape (N,)
        Wavelength vector for the selected range
    WLight : dict
        Full structure loaded from the .mat file
    idx_min : int
        Index of wl_min in the full lambda vector
    idx_max : int
        Index of wl_max in the full lambda vector
    lambda_full : ndarray
        Full wavelength vector (not truncated)
    """
    WLight = loadmat(path, squeeze_me=True, struct_as_record=False)

    # Full wavelength vector
    if 'lambda' in WLight:
        lambda_full = np.asarray(WLight['lambda']).flatten().astype(float)
    elif 'raw_lambda' in WLight:
        lambda_full = np.asarray(WLight['raw_lambda']).flatten().astype(float)
    else:
        raise KeyError(f"No wavelength vector in {path}")

    # Indices of the range of interest
    idx_min = _find_wl_index(lambda_full, wl_min)
    idx_max = _find_wl_index(lambda_full, wl_max)
    lambda_wl = lambda_full[idx_min:idx_max + 1]

    # Acquisition metadata
    donnees = WLight['donnees_acq_brut']
    acq_data = np.asarray(donnees.data if hasattr(donnees, 'data') else donnees)

    signal_brut = np.asarray(WLight['signal_brut'])

    # WLallummee = columns 4,5,6 MATLAB (3,4,5 Python)
    # idxShort = col 0 of WLallummee = col 3 of acq_data
    # idxLong  = col 1 of WLallummee = col 4 of acq_data
    led_cols = acq_data[:, 3:6]
    idxShort = np.where(led_cols[:, 0] == 1)[0]
    idxLong = np.where(led_cols[:, 1] == 1)[0]

    def _extract(idx_array):
        """Background subtraction for each acquisition."""
        spectra = []
        for i in idx_array:
            if i > 0:
                diff = signal_brut[i, idx_min:idx_max + 1] - signal_brut[i - 1, idx_min:idx_max + 1]
                spectra.append(diff)
        if not spectra:
            return np.zeros((1, len(lambda_wl)))
        return np.array(spectra, dtype=float)

    allShort = _extract(idxShort)
    allLong = _extract(idxLong)

    if return_all:
        WLShort, WLLong = allShort, allLong
    else:
        WLShort, WLLong = np.mean(allShort, axis=0), np.mean(allLong, axis=0)

    return WLShort, WLLong, lambda_wl, WLight, idx_min, idx_max, lambda_full
