"""
Python translation of corrected_fluo_LS_WL_1.m (Antoine Uzel, 24/08/2024)

Complete fluorescence correction pipeline with optional extraction of
optical properties by fitting the Kim model.

Note on optical correction:
  - In the original MATLAB code, R_short_extracted is forced to 1 (line 75)
    even after the reflectance fit. The fit is therefore done but not applied.
  - In this Python version, the behavior is controlled by `apply_optical_correction`:
    * False (default) -> R = 1, identical to the current MATLAB behavior
    * True  -> uses R_short_extracted from the fit ("full" behavior)
"""

import os
import numpy as np
from scipy.io import loadmat
from scipy.optimize import least_squares
from scipy.interpolate import interp1d

from .load_fluo import load_fluo
from .load_reflectance import load_reflectance
from .calibration_spectralon import calibration_spectralon
from .reflectance_kim import reflectance_kim
from .fluo_model import fluo_exp_385, extract_components


def _search_file(start_dir, filename, max_levels=4):
    """Searches for a file in start_dir and its parents up to max_levels."""
    current = start_dir
    for _ in range(max_levels):
        candidate = os.path.join(current, filename)
        if os.path.isfile(candidate):
            return candidate
        parent = os.path.dirname(current)
        if parent == current:
            break
        current = parent
    return None


def _load_absorption_spectra(search_dir, lambda_wl):
    """
    Loads the biological absorption spectra and interpolates them onto lambda_wl.

    Returns ua of shape (5, N):
      ua[0] : HbO2
      ua[1] : Hb
      ua[2] : mua_fat (fat)
      ua[3] : oxidized cytochrome b
      ua[4] : reduced cytochrome b
    """
    N = len(lambda_wl)
    ua = np.zeros((5, N))

    CHb = 150.0    # Hb concentration, g/L
    MHb = 64500.0  # Hb molar mass, g/mol

    # HbO2 and Hb
    hb_file = _search_file(search_dir, 'HbO2 Hb visible.txt')
    if hb_file:
        try:
            hb_data = np.loadtxt(hb_file)
            interp_HbO2 = interp1d(hb_data[:, 0], 2.303 * hb_data[:, 1] / MHb,
                                    bounds_error=False, fill_value=0.0)
            interp_Hb = interp1d(hb_data[:, 0], 2.303 * hb_data[:, 2] / MHb,
                                  bounds_error=False, fill_value=0.0)
            ua[0, :] = interp_HbO2(lambda_wl)
            ua[1, :] = interp_Hb(lambda_wl)
        except Exception as e:
            print(f"[Warning] Could not load HbO2/Hb: {e}")

    # Fat (mua_fat) — the lambda field is a Python keyword, use getattr
    fat_file = _search_file(search_dir, 'mua_fat.mat')
    if fat_file:
        try:
            fat_data = loadmat(fat_file, squeeze_me=True, struct_as_record=False)
            absorb = fat_data['absorbeur']
            # Try the possible names for the wavelength field
            fat_lam = None
            for field in ['lambda', 'wl', 'wavelength', 'lambda_wl']:
                val = getattr(absorb, field, None)
                if val is not None:
                    fat_lam = np.asarray(val).flatten()
                    break
            fat_mua = np.asarray(absorb.mua).flatten()
            if fat_lam is not None:
                interp_fat = interp1d(fat_lam, fat_mua, bounds_error=False, fill_value=0.0)
                ua[2, :] = interp_fat(lambda_wl)
        except Exception as e:
            print(f"[Warning] Could not load mua_fat: {e}")

    # Oxidized cytochrome b
    cyt_oxi_file = _search_file(search_dir, 'cyt b oxidized.txt')
    if cyt_oxi_file:
        try:
            cyt_data = np.loadtxt(cyt_oxi_file)
            interp_cyt = interp1d(cyt_data[:, 0], cyt_data[:, 1],
                                   bounds_error=False, fill_value=0.0)
            ua[3, :] = interp_cyt(lambda_wl)
        except Exception as e:
            print(f"[Warning] Could not load oxidized cyt b: {e}")

    # Reduced cytochrome b
    cyt_red_file = _search_file(search_dir, 'cyt b reduced.txt')
    if cyt_red_file:
        try:
            cyt_data = np.loadtxt(cyt_red_file)
            interp_cyt = interp1d(cyt_data[:, 0], cyt_data[:, 1],
                                   bounds_error=False, fill_value=0.0)
            ua[4, :] = interp_cyt(lambda_wl)
        except Exception as e:
            print(f"[Warning] Could not load reduced cyt b: {e}")

    return ua


def _fit_reflectance(R_measured, lambda_wl, ua, n_probe, rho=240e-4, n_medium=1.4):
    """
    Fits the Kim model to the measured reflectance.

    Parameters x (8 elements):
      x[0]=c1, x[1]=c2, x[2]=d, x[3]=c_blood, x[4]=StO2, x[5]=unused, x[6]=cyt_oxi, x[7]=cyt_red

    Returns (x_opt, R_fitted).
    """
    x0 = np.array([10.0, 1.0, 0.0, 0.9, 0.5, 0.5, 1.0, 1.0])
    lb = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    ub = np.array([100.0, 4.0, 1.0, 70.0, 1.0, 100.0, np.inf, np.inf])

    def residual(x):
        pred = reflectance_kim(x, rho, ua, lambda_wl, n_medium, n_probe)
        return pred - R_measured

    try:
        result = least_squares(residual, x0, bounds=(lb, ub), method='trf',
                               max_nfev=int(1e6), ftol=1e-10, xtol=1e-10, gtol=1e-10)
        x_opt = result.x
        R_fitted = reflectance_kim(x_opt, rho, ua, lambda_wl, n_medium, n_probe)
        return x_opt, R_fitted, result.cost
    except Exception as e:
        print(f"[Warning] Reflectance fit failed: {e}")
        return x0, np.ones_like(R_measured), np.inf


def _fit_fluorescence(spectrum, lambda_wl, fluorophores, R_correction,
                      active_fluorophores=None, ref_spectra=None):
    """
    Fits the fluorescence model to the measured spectrum.

    Full parameters (15 elements):
      p[0]  : alpha
      p[1-6]: amplitudes FAD, NADH, FMN, Lipo, PpIX_636, PpIX_620
      p[7-14]: Gaussian parameters (means/std devs FMN, Lipo, PpIX_636, PpIX_620)
               — fixed automatically if a reference spectrum is provided in ref_spectra

    ref_spectra : optional dict {'FMN': array, 'Lipo': array, 'PpIX_636': array, 'PpIX_620': array}
        When a key is present, the spectrum replaces the Gaussian and its shape
        parameters are fixed (only the amplitude remains free).
    """
    if active_fluorophores is None:
        active_fluorophores = {k: True for k in
                               ['FAD', 'NADH', 'FMN', 'Lipo', 'PpIX_636', 'PpIX_620']}
    if ref_spectra is None:
        ref_spectra = {}

    lb = np.array([0.,   0.,   0.,   0.,    0.,    0.,    0.,   494., 14., 589.,  9.,  636., 5.5, 618., 7.5], dtype=float)
    ub = np.array([1., 1000., 1000., 100., 1000., 1000., 1000., 496., 16., 591., 11.,  638., 7.5, 620.,   9.], dtype=float)
    x0 = np.array([0.1, 10.,  10.,   1.,    5.,    1.,    1.,   495., 15., 590., 10.,  637., 6.5, 619., 8.25], dtype=float)

    gaussian_params = {
        'FMN':      (3, 7,  8),
        'Lipo':     (4, 9,  10),
        'PpIX_636': (5, 11, 12),
        'PpIX_620': (6, 13, 14),
    }

    for name, amp_idx in [('FAD', 1), ('NADH', 2)]:
        if not active_fluorophores.get(name, True):
            lb[amp_idx] = ub[amp_idx] = x0[amp_idx] = 0.

    for name, (amp_idx, mean_idx, std_idx) in gaussian_params.items():
        if not active_fluorophores.get(name, True):
            # Fluorophore disabled: amplitude and Gaussian parameters fixed
            lb[amp_idx]  = ub[amp_idx]  = x0[amp_idx]  = 0.
            lb[mean_idx] = ub[mean_idx] = x0[mean_idx] = (lb[mean_idx] + ub[mean_idx]) / 2.
            lb[std_idx]  = ub[std_idx]  = x0[std_idx]  = (lb[std_idx]  + ub[std_idx])  / 2.
        elif name in ref_spectra:
            # Reference spectrum provided: Gaussian parameters unnecessary -> fixed
            lb[mean_idx] = ub[mean_idx] = x0[mean_idx] = (lb[mean_idx] + ub[mean_idx]) / 2.
            lb[std_idx]  = ub[std_idx]  = x0[std_idx]  = (lb[std_idx]  + ub[std_idx])  / 2.

    free = lb < ub
    fixed_values = x0.copy()
    x0_free = x0[free]
    lb_free = lb[free]
    ub_free = ub[free]

    R_is_scalar = np.isscalar(R_correction) or (
        isinstance(R_correction, np.ndarray) and R_correction.size == 1)

    def _build_full(p_free):
        p = fixed_values.copy()
        p[free] = p_free
        return p

    def residual(p_free):
        p = _build_full(p_free)
        alpha = p[0]
        fluo_params = p[1:]
        spectrum_fit = fluo_exp_385(fluo_params, fluorophores, lambda_wl, ref_spectra)
        correction = (float(R_correction) ** alpha) if R_is_scalar else (R_correction ** alpha)
        return spectrum_fit * correction - spectrum

    try:
        result = least_squares(residual, x0_free, bounds=(lb_free, ub_free), method='trf',
                               max_nfev=int(1e8), ftol=1e-6, xtol=1e-10)
        return _build_full(result.x), result.cost
    except Exception as e:
        print(f"[Warning] Fluorescence fit failed: {e}")
        return fixed_values, np.inf


def corrected_fluo_ls_wl_1(
    path,
    file,
    name_spectralon='spectralon_000_diffuseReflectances',
    is_small_spectralon=True,
    min_wl_reflectance=480,
    max_wl_reflectance=645,
    min_wl_fluo=480,
    max_wl_fluo=645,
    fluorophores_385=None,
    fluorophores_405=None,
    apply_optical_correction=False,
    active_fluorophores=None,
    spectralon_theo_path=None,
    scale_path=None,
    ref_spectra_385=None,
    ref_spectra_405=None,
):
    """
    Complete pipeline: loads the fluorescence, optionally corrects it using the
    optical properties, fits the fluorescence model, and decomposes it into components.

    Translated from corrected_fluo_LS_WL_1.m (Antoine Uzel, 24/08/2024).

    Parameters
    ----------
    path : str
        Directory containing the data files
    file : str
        Fluorescence file name (*fluo.mat)
    name_spectralon : str
        Base name of the spectralon file (without extension)
    is_small_spectralon : bool
        True -> x0.8 correction for the small spectralon
    min_wl_reflectance, max_wl_reflectance : float
        Wavelength range for reflectance (nm)
    min_wl_fluo, max_wl_fluo : float
        Wavelength range for fluorescence (nm)
    fluorophores_385 : ndarray, shape (2, N)
        FAD and NADH reference spectra for 385 nm excitation
    fluorophores_405 : ndarray, shape (2, N)
        FAD and NADH reference spectra for 405 nm excitation
    apply_optical_correction : bool
        False -> R = 1 (current MATLAB behavior)
        True  -> uses R extracted from the reflectance fit
    active_fluorophores : dict, optional
        {'FAD': True/False, 'NADH': True/False, ...}
    spectralon_theo_path : str, optional
        Path to Reflectance_values_array.txt
    scale_path : str, optional
        Path to scale_new_theory_intralipids.mat

    Returns
    -------
    result : dict with keys:
        'S385total', 'S405total'         : raw fluorescence (N,)
        'S385_corrected', 'S405_corrected': corrected fluorescence (N,)
        'fluorophore'                     : dict of decomposed components
        'lambda_fluo'                     : fluorescence wavelengths (N,)
        'lambda_reflectance'              : reflectance wavelengths (N,) or None
        'res_385', 'res_405'              : fit parameters (15,)
        'R_short'                         : reflectance used for correction
        'optical_params'                  : Kim parameters (8,) if optical correction applied
    """
    if active_fluorophores is None:
        active_fluorophores = {k: True for k in
                               ['FAD', 'NADH', 'FMN', 'Lipo', 'PpIX_636', 'PpIX_620']}

    result = {}

    # --- 1. Load fluorescence ---
    fluo_path = os.path.join(path, file)
    S385total, S405total, lambda_fluo = load_fluo(fluo_path, min_wl_fluo, max_wl_fluo)
    result['S385total'] = S385total
    result['S405total'] = S405total
    result['lambda_fluo'] = lambda_fluo

    R_short_extracted = 1.0  # default value (as in MATLAB)
    result['lambda_reflectance'] = None
    result['optical_params'] = None

    if apply_optical_correction:
        # --- 2. Load reflectance ---
        # strrep(file, 'fluo', 'diffuseReflectances') as in MATLAB
        file_wl = file.replace('fluo', 'diffuseReflectances')
        reflectance_path = os.path.join(path, file_wl)

        if not os.path.isfile(reflectance_path):
            raise FileNotFoundError(
                f"Diffuse reflectance file not found: {reflectance_path}\n"
                "Disable optical correction or check the file."
            )

        WLShort, WLLong, lambda_wl, WLight, idx_min, idx_max, lambda_full = \
            load_reflectance(reflectance_path, min_wl_reflectance, max_wl_reflectance)
        result['lambda_reflectance'] = lambda_wl

        # --- 3. Calibration scale factors ---
        if scale_path is None:
            scale_path = _search_file(path, 'scale_new_theory_intralipids.mat')
        if scale_path is None:
            raise FileNotFoundError("scale_new_theory_intralipids.mat not found")

        scale_data = loadmat(scale_path, squeeze_me=True, struct_as_record=False)
        Scale = scale_data['Scale']
        WL_Short = np.asarray(Scale.Short.WL).flatten()
        Sy_Short = np.asarray(Scale.Short.Sy).flatten()
        n_Short_arr = np.asarray(Scale.Short.n).flatten()

        scale_short = interp1d(WL_Short, Sy_Short, bounds_error=False,
                               fill_value='extrapolate')(lambda_wl)
        # n_probe interpolated over lambda_full then truncated as in MATLAB
        n_short_full = interp1d(WL_Short, n_Short_arr, bounds_error=False,
                                fill_value='extrapolate')(lambda_full)
        n_probe_wl = n_short_full[idx_min:idx_max + 1]

        # --- 4. Spectralon calibration ---
        # First check whether the .mat file contains a 'name_spectralon' field
        if 'name_spectralon' in WLight:
            spectralon_exp_path = os.path.join(path, str(WLight['name_spectralon']))
        else:
            spectralon_exp_path = os.path.join(path, name_spectralon + '.mat')

        if spectralon_theo_path is None:
            spectralon_theo_path = _search_file(path, 'Reflectance_values_array.txt')
        if spectralon_theo_path is None:
            raise FileNotFoundError("Reflectance_values_array.txt not found")

        SsourceShort, _, lambda_spec = calibration_spectralon(
            spectralon_exp_path, spectralon_theo_path, is_small_spectralon
        )
        # Truncate to the reflectance range
        SsourceShort_wl = SsourceShort[idx_min:idx_max + 1]

        # --- 5. Measured reflectance ---
        R_Short = WLShort / (SsourceShort_wl * scale_short)

        # --- 6. Absorption spectra ---
        ua = _load_absorption_spectra(path, lambda_wl)

        # --- 7. Kim fit ---
        rho_short = 240e-4  # 2.4 mm in meters
        n_medium = 1.4
        x_opt, R_fitted, _ = _fit_reflectance(
            R_Short, lambda_wl, ua, n_probe_wl, rho=rho_short, n_medium=n_medium
        )
        result['optical_params'] = x_opt

        # R_Short (measured) is used for the fluorescence correction.
        # The Kim model is used only to extract the optical properties (StO2, c_blood).
        # (Equivalent to line 74 of the original MATLAB, before the R=1 override.)
        R_short_extracted = R_Short
        result['R_short'] = R_short_extracted
        result['R_short_fitted'] = R_fitted
    else:
        result['R_short'] = R_short_extracted

    # R must be on the same grid as the fluorescence before being passed to the fit.
    # If R_short_extracted is an array (reflectance grid), interpolate it onto lambda_fluo.
    R_for_fit = R_short_extracted
    if not np.isscalar(R_short_extracted):
        _interp_R = interp1d(result['lambda_reflectance'], R_short_extracted,
                             bounds_error=False, fill_value='extrapolate')
        R_for_fit = _interp_R(lambda_fluo)

    # --- 8. Interpolate reference spectra onto lambda_fluo ---
    # The provided spectra may have a slightly different grid (different file).
    # They are re-interpolated onto lambda_fluo to guarantee dimensional compatibility.
    # If the reference spectra carry their own lambda axis, it must be passed.
    # By convention we accept either a 1D array (emission only, unknown grid ->
    # assumed to span the same range as lambda_fluo), or a tuple (lambda, values).
    def _interp_ref_smart(ref_dict, lam_target):
        if not ref_dict:
            return {}
        out = {}
        for name, val in ref_dict.items():
            if isinstance(val, (list, tuple)) and len(val) == 2:
                lam_src, arr = np.asarray(val[0], dtype=float), np.asarray(val[1], dtype=float)
            else:
                arr = np.asarray(val, dtype=float).flatten()
                lam_src = np.linspace(lam_target[0], lam_target[-1], len(arr))
            out[name] = interp1d(lam_src, arr,
                                 bounds_error=False, fill_value=0.0)(lam_target)
        return out

    ref385 = _interp_ref_smart(ref_spectra_385, lambda_fluo)
    ref405 = _interp_ref_smart(ref_spectra_405, lambda_fluo)

    # --- 9. Fluorescence fit ---
    if fluorophores_385 is None or fluorophores_405 is None:
        raise ValueError("The reference spectra fluorophores_385 and fluorophores_405 "
                         "must be provided.")

    res_385, cost_385 = _fit_fluorescence(
        S385total, lambda_fluo, fluorophores_385, R_for_fit, active_fluorophores, ref385)
    res_405, cost_405 = _fit_fluorescence(
        S405total, lambda_fluo, fluorophores_405, R_for_fit, active_fluorophores, ref405)

    result['res_385'] = res_385
    result['res_405'] = res_405

    # --- 9. Fluorescence correction ---
    alpha_385 = res_385[0]
    alpha_405 = res_405[0]

    # R_for_fit is already interpolated onto lambda_fluo
    if np.isscalar(R_for_fit):
        R_val = float(R_for_fit)
        S385_corrected = S385total / (R_val ** alpha_385) if R_val != 0 else S385total
        S405_corrected = S405total / (R_val ** alpha_405) if R_val != 0 else S405total
    else:
        R_safe = np.maximum(R_for_fit, 1e-10)
        S385_corrected = S385total / (R_safe ** alpha_385)
        S405_corrected = S405total / (R_safe ** alpha_405)

    result['S385_corrected'] = S385_corrected
    result['S405_corrected'] = S405_corrected

    # --- 10. Decomposition into components ---
    fluo_params_385 = res_385[1:]  # without alpha
    fluo_params_405 = res_405[1:]

    comp_385 = extract_components(fluo_params_385, fluorophores_385, lambda_fluo, ref385)
    comp_405 = extract_components(fluo_params_405, fluorophores_405, lambda_fluo, ref405)

    fluorophore = {}
    name_map = {'FAD': 'flavine', 'NADH': 'NADH', 'FMN': 'gaussian',
                'Lipo': 'lipo', 'PpIX_636': 'PpIX_636', 'PpIX_620': 'PpIX_620'}
    for key in comp_385:
        fluorophore[f'{name_map.get(key, key)}_385_exp'] = comp_385[key]
        fluorophore[f'{name_map.get(key, key)}_405_exp'] = comp_405[key]

    # Ratios redox
    FAD_385_tot = np.sum(comp_385['FAD'])
    NADH_385_tot = np.sum(comp_385['NADH'])
    FAD_405_tot = np.sum(comp_405['FAD'])
    NADH_405_tot = np.sum(comp_405['NADH'])

    fluorophore['FAD_385_tot'] = FAD_385_tot
    fluorophore['NADH_385_tot'] = NADH_385_tot
    fluorophore['FAD_405_tot'] = FAD_405_tot
    fluorophore['NADH_405_tot'] = NADH_405_tot
    fluorophore['redox_385'] = FAD_385_tot / (FAD_385_tot + NADH_385_tot) \
        if (FAD_385_tot + NADH_385_tot) > 0 else 0.0
    fluorophore['redox_405'] = FAD_405_tot / (FAD_405_tot + NADH_405_tot) \
        if (FAD_405_tot + NADH_405_tot) > 0 else 0.0
    fluorophore['StO2'] = float(result['optical_params'][4]) \
        if result['optical_params'] is not None else None
    fluorophore['c_blood'] = float(result['optical_params'][3]) \
        if result['optical_params'] is not None else None

    result['fluorophore'] = fluorophore

    return result
