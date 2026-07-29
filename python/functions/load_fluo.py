"""
Python translation of load_fluo.m
Loads and processes fluorescence data from a *fluo.mat file.

In the acquisition data:
  - donnees_acq_brut.data(:, 4:6) [MATLAB 1-indexed] -> columns 3,4,5 [Python 0-indexed]
  - Column 3 (WLallummee col 0) = 385 nm laser
  - Column 4 (WLallummee col 1) = 405 nm laser
  - Column 6 = Power
  - Column 7 = Acquisition time (ms)
"""

import numpy as np
from scipy.io import loadmat


def _find_wl_index(lambda_arr, wl_target):
    """Returns the index closest to wl_target in lambda_arr."""
    return int(np.argmin((lambda_arr - wl_target) ** 2))


def _load_mat_safe(path):
    """Loads a .mat file with handling of nested structures."""
    return loadmat(path, squeeze_me=True, struct_as_record=False)


def load_fluo(path, wl_min=480, wl_max=645, return_all=False):
    """
    Loads fluorescence data from a *fluo.mat file.

    Translated from load_fluo.m (Antoine Uzel).

    Parameters
    ----------
    path : str
        Full path to the *fluo.mat file
    wl_min : float
        Minimum wavelength of interest (nm)
    wl_max : float
        Maximum wavelength of interest (nm)
    return_all : bool
        If True, returns all acquisitions (2D) instead of the average (1D)

    Returns
    -------
    fluo_385 : ndarray
        Fluorescence excitation 385 nm — (N,) if return_all=False, (M, N) otherwise
    fluo_405 : ndarray
        Fluorescence excitation 405 nm — (N,) if return_all=False, (M, N) otherwise
    lambda_wl : ndarray
        Wavelength vector (nm), shape (N,)
    """
    data = _load_mat_safe(path)

    # Wavelength vector
    if 'lambda' in data:
        lambda_full = np.asarray(data['lambda']).flatten().astype(float)
    elif 'raw_lambda' in data:
        lambda_full = np.asarray(data['raw_lambda']).flatten().astype(float)
    else:
        raise KeyError(f"No wavelength vector found in {path}")

    idx_min = _find_wl_index(lambda_full, wl_min)
    idx_max = _find_wl_index(lambda_full, wl_max)
    lambda_wl = lambda_full[idx_min:idx_max + 1]

    # Acquisition metadata
    donnees = data['donnees_acq_brut']
    acq_data = np.asarray(donnees.data if hasattr(donnees, 'data') else donnees)

    signal_brut = np.asarray(data['signal_brut'])

    # LEDallummee = columns 4,5,6 of data (MATLAB 1-indexed) -> 3,4,5 (Python 0-indexed)
    # idx385: column 1 of LEDallummee = column 3 of data
    # idx405: column 2 of LEDallummee = column 4 of data
    led_cols = acq_data[:, 3:6]
    idx385 = np.where(led_cols[:, 0] == 1)[0]
    idx405 = np.where(led_cols[:, 1] == 1)[0]

    if len(idx385) == 0:
        raise ValueError("No acquisition with the 385 nm laser found in this file.")
    if len(idx405) == 0:
        raise ValueError("No acquisition with the 405 nm laser found in this file.")

    # Power and acquisition time (normalization)
    power385 = float(acq_data[idx385[0], 6])
    power405 = float(acq_data[idx405[0], 6])
    time_acq = float(acq_data[0, 7]) / 1000.0  # ms -> seconds

    norm385 = power385 * time_acq if power385 * time_acq != 0 else 1.0
    norm405 = power405 * time_acq if power405 * time_acq != 0 else 1.0

    def _extract(idx_array, norm):
        """Background subtraction + normalization for each acquisition."""
        spectra = []
        for i in idx_array:
            if i > 0:
                diff = signal_brut[i, :] - signal_brut[i - 1, :]
                spectra.append(diff[idx_min:idx_max + 1] / norm)
        if not spectra:
            return np.zeros((1, len(lambda_wl)))
        return np.array(spectra, dtype=float)

    all385 = _extract(idx385, norm385)
    all405 = _extract(idx405, norm405)

    if return_all:
        return all385, all405, lambda_wl
    else:
        return np.mean(all385, axis=0), np.mean(all405, axis=0), lambda_wl
