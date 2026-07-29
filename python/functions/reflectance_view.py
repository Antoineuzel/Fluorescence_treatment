"""
Loading of diffuse reflectance for visualization (GUI tab + scripts).

Deliberately separates the visualization (this module) from the optical
correction using the Kim model (functions/corrected_fluo.py): here no fit is
performed, we only compute the raw signal and, if the calibration files are
available, the reflectance calibrated by the spectralon.
"""

import os
import numpy as np
from scipy.io import loadmat
from scipy.interpolate import interp1d

from .load_reflectance import load_reflectance
from .calibration_spectralon import calibration_spectralon


def load_reflectance_view(fluo_path, spectralon_theo_path=None, scale_path=None,
                           is_small_spectralon=False, wl_min=480, wl_max=645,
                           name_spectralon='spectralon_000_diffuseReflectances',
                           return_all=False):
    """
    Loads the short/long wavelength reflectance associated with a *fluo.mat
    file, for display only (no fit).

    Parameters
    ----------
    fluo_path : str
        Full path to the currently selected *fluo.mat file.
    spectralon_theo_path : str, optional
        Path to Reflectance_values_array.txt. If absent, only the raw
        signal is returned (no calibration).
    scale_path : str, optional
        Path to scale_new_theory_intralipids.mat.
    is_small_spectralon : bool
        If True, corrects for the small spectralon (factor 0.8).
    wl_min, wl_max : float
        Wavelength range (nm).
    name_spectralon : str
        Base name of the experimental spectralon file expected in the same
        folder as fluo_path, if the reflectance .mat file does not itself
        reference its spectralon file.
    return_all : bool
        If True, WLShort/WLLong contain all acquisitions (2D) instead of
        their average (1D). The calibrated reflectance is still computed on
        the average (calibration_spectralon only works on averages).

    Returns
    -------
    dict with keys:
        'reflectance_path' : str — path to the *diffuseReflectances.mat file used
        'lambda'           : ndarray (N,)
        'WLShort', 'WLLong'          : ndarray (N,) or (M, N) — raw signal
        'R_Short', 'R_Long'          : ndarray (N,) or None — calibrated reflectance (average)
        'calibration_error'          : str or None — reason for calibration failure
    """
    directory = os.path.dirname(fluo_path)
    fluo_file = os.path.basename(fluo_path)
    reflectance_file = fluo_file.replace('fluo', 'diffuseReflectances')
    reflectance_path = os.path.join(directory, reflectance_file)

    if not os.path.isfile(reflectance_path):
        raise FileNotFoundError(
            f"Reflectance file not found: {reflectance_file}\n"
            f"(expected in {directory})"
        )

    WLShort, WLLong, lambda_wl, WLight, idx_min, idx_max, lambda_full = \
        load_reflectance(reflectance_path, wl_min, wl_max, return_all=return_all)

    result = {
        'reflectance_path': reflectance_path,
        'lambda': lambda_wl,
        'WLShort': WLShort,
        'WLLong': WLLong,
        'R_Short': None,
        'R_Long': None,
        'calibration_error': None,
    }

    if spectralon_theo_path is None or not os.path.isfile(spectralon_theo_path):
        result['calibration_error'] = "Reflectance_values_array.txt not found."
        return result
    if scale_path is None or not os.path.isfile(scale_path):
        result['calibration_error'] = "scale_new_theory_intralipids.mat not found."
        return result

    try:
        if 'name_spectralon' in WLight:
            spectralon_exp_path = os.path.join(directory, str(WLight['name_spectralon']))
        else:
            spectralon_exp_path = os.path.join(directory, name_spectralon + '.mat')

        if not os.path.isfile(spectralon_exp_path):
            result['calibration_error'] = (
                f"Spectralon file not found: {os.path.basename(spectralon_exp_path)}"
            )
            return result

        SsourceShort, SsourceLong, _ = calibration_spectralon(
            spectralon_exp_path, spectralon_theo_path, is_small_spectralon
        )
        SsourceShort_wl = SsourceShort[idx_min:idx_max + 1]
        SsourceLong_wl = SsourceLong[idx_min:idx_max + 1]

        scale_data = loadmat(scale_path, squeeze_me=True, struct_as_record=False)
        Scale = scale_data['Scale']

        scale_short = interp1d(np.asarray(Scale.Short.WL).flatten(),
                                np.asarray(Scale.Short.Sy).flatten(),
                                bounds_error=False, fill_value='extrapolate')(lambda_wl)
        scale_long = interp1d(np.asarray(Scale.Long.WL).flatten(),
                               np.asarray(Scale.Long.Sy).flatten(),
                               bounds_error=False, fill_value='extrapolate')(lambda_wl)

        result['R_Short'] = WLShort / (SsourceShort_wl * scale_short)
        result['R_Long'] = WLLong / (SsourceLong_wl * scale_long)
    except Exception as e:
        result['calibration_error'] = f"Calibration failed: {e}"

    return result
