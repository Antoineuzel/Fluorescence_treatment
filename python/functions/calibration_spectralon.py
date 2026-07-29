"""
Python translation of calibration_spectralon.m and calibration_spectralon_from_data.m

Computes the calibrated source spectrum (SsourceShort, SsourceLong) from
experimental spectralon measurements and theoretical reference values.

Formula: SsourceShort = sum(SpecexpShort) / (theoreticalSpectralon * nb_acquisitions)
"""

import numpy as np
from scipy.io import loadmat
from scipy.interpolate import interp1d


def _load_mat_safe(path):
    return loadmat(path, squeeze_me=True, struct_as_record=False)


def _load_reflectance_values(filepath):
    """
    Loads the Reflectance_values_array.txt file.
    Format: 2 header lines, then 2 columns (wavelength, reflectance).
    """
    try:
        data = np.loadtxt(filepath, skiprows=2)
    except Exception:
        data = np.loadtxt(filepath, skiprows=1)
    return data


def calibration_spectralon(file_spectralon_exp, file_spectralon_theo, is_small_spectralon=False):
    """
    Calibrates the source spectrum from a spectralon measurement.

    Translated from calibration_spectralon.m (uses signal_brut).

    Parameters
    ----------
    file_spectralon_exp : str
        Path to the experimental spectralon .mat file
    file_spectralon_theo : str
        Path to the theoretical reflectance text file (Reflectance_values_array.txt)
    is_small_spectralon : bool
        If True, divides the source spectrum by 0.8 (small spectralon correction)

    Returns
    -------
    SsourceShort : ndarray
        Calibrated source spectrum — short wavelength (vector over the full lambda range)
    SsourceLong : ndarray
        Calibrated source spectrum — long wavelength
    lambda_exp : ndarray
        Wavelength vector of the experimental spectralon
    """
    spec_exp = _load_mat_safe(file_spectralon_exp)

    # Wavelength vector
    if 'lambda' in spec_exp:
        lambda_exp = np.asarray(spec_exp['lambda']).flatten().astype(float)
    elif 'raw_lambda' in spec_exp:
        lambda_exp = np.asarray(spec_exp['raw_lambda']).flatten().astype(float)
    else:
        raise KeyError(f"No lambda vector in {file_spectralon_exp}")

    donnees = spec_exp['donnees_acq_brut']
    acq_data = np.asarray(donnees.data if hasattr(donnees, 'data') else donnees)
    signal_brut = np.asarray(spec_exp['signal_brut'])

    # LED columns: 3,4,5 (Python) = 4,5,6 (MATLAB)
    led_cols = acq_data[:, 3:6]
    idxShort = np.where(led_cols[:, 0] == 1)[0]
    idxLong = np.where(led_cols[:, 1] == 1)[0]

    def _subtract_background(idx_array):
        spectra = []
        for i in idx_array:
            if i > 0:
                spectra.append(signal_brut[i, :] - signal_brut[i - 1, :])
        return np.array(spectra, dtype=float) if spectra else np.zeros((1, signal_brut.shape[1]))

    SpecexpShort = _subtract_background(idxShort)
    SpecexpLong = _subtract_background(idxLong)

    # Theoretical reflectance of the spectralon
    theo_data = _load_reflectance_values(file_spectralon_theo)
    theo_wl = theo_data[:, 0]
    theo_refl = theo_data[:, 1]

    interp_func = interp1d(theo_wl, theo_refl, bounds_error=False, fill_value='extrapolate')
    spec_theo = interp_func(lambda_exp)

    # Source spectrum: sum of acquisitions normalized by the theoretical reflectance
    n_short = max(len(idxShort), 1)
    n_long = max(len(idxLong), 1)

    SsourceShort = np.sum(SpecexpShort, axis=0) / (spec_theo * n_short)
    SsourceLong = np.sum(SpecexpLong, axis=0) / (spec_theo * n_long)

    # Small spectralon correction: divide the source spectrum by 0.8
    # (the small spectralon reflects 80% of what the large spectralon reflects)
    if is_small_spectralon:
        SsourceShort = SsourceShort / 0.8
        SsourceLong = SsourceLong / 0.8

    return SsourceShort, SsourceLong, lambda_exp


def calibration_spectralon_from_data(file_wl, file_spectralon_theo, is_small_spectralon=False):
    """
    Alternative calibration using the 'spectralon' field pre-extracted from the .mat file.

    Translated from calibration_spectralon_from_data.m.
    Used when the spectralon data is stored in a separate field.
    """
    spec_exp = _load_mat_safe(file_wl)

    if 'lambda' in spec_exp:
        lambda_exp = np.asarray(spec_exp['lambda']).flatten().astype(float)
    else:
        lambda_exp = np.asarray(spec_exp['raw_lambda']).flatten().astype(float)

    donnees = spec_exp['donnees_acq_brut']
    acq_data = np.asarray(donnees.data if hasattr(donnees, 'data') else donnees)

    # Uses the 'spectralon' field instead of 'signal_brut'
    spectralon_data = np.asarray(spec_exp['spectralon'])

    led_cols = acq_data[:, 3:6]
    idxShort = np.where(led_cols[:, 0] == 1)[0]
    idxLong = np.where(led_cols[:, 1] == 1)[0]

    def _subtract_background(idx_array):
        spectra = []
        for i in idx_array:
            if i > 0:
                spectra.append(spectralon_data[i, :] - spectralon_data[i - 1, :])
        return np.array(spectra, dtype=float) if spectra else np.zeros((1, spectralon_data.shape[1]))

    SpecexpShort = _subtract_background(idxShort)
    SpecexpLong = _subtract_background(idxLong)

    theo_data = _load_reflectance_values(file_spectralon_theo)
    interp_func = interp1d(theo_data[:, 0], theo_data[:, 1], bounds_error=False, fill_value='extrapolate')
    spec_theo = interp_func(lambda_exp)

    n_short = max(len(idxShort), 1)
    n_long = max(len(idxLong), 1)

    SsourceShort = np.sum(SpecexpShort, axis=0) / (spec_theo * n_short)
    SsourceLong = np.sum(SpecexpLong, axis=0) / (spec_theo * n_long)

    if is_small_spectralon:
        SsourceShort = SsourceShort / 0.8
        SsourceLong = SsourceLong / 0.8

    return SsourceShort, SsourceLong, lambda_exp
