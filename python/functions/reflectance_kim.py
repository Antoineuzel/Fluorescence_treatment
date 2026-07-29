"""
Python translation of reflectance_Kim.m (Arthur / Antoine Uzel)

Radiative transfer model (diffusion approximation) of Kim et al.
Computes the diffuse reflectance from the optical properties of the tissue.

Parameters of vector x (8 elements, Python 0-indexed):
  x[0] : c1 — scattering amplitude
  x[1] : c2 — scattering exponent (power law)
  x[2] : d  — Mie / Rayleigh fraction
  x[3] : c_blood — blood concentration (g/L)
  x[4] : StO2   — oxygen saturation (0-1)
  x[5] : (unused, mua_fat coefficient commented out in MATLAB)
  x[6] : oxidized cytochrome b coefficient
  x[7] : reduced cytochrome b coefficient
"""

import numpy as np


def reflectance_kim(x, rho, ua, lambda_wl, n_medium, n_probe):
    """
    Kim diffuse reflectance model (2-source diffusion approximation).

    Translated from reflectance_Kim.m.

    Parameters
    ----------
    x : array-like, length 8
        Model parameters (see module header)
    rho : float
        Source-detector distance (in meters), typically 2.4e-3 m
    ua : ndarray, shape (5, N)
        Absorption coefficients of the chromophores:
          - ua[0] : HbO2
          - ua[1] : deoxygenated Hb
          - ua[2] : fat (mua_fat) — unused in the final expression
          - ua[3] : oxidized cytochrome b
          - ua[4] : reduced cytochrome b
    lambda_wl : ndarray, shape (N,)
        Wavelength vector (nm) — pre-truncated to the range of interest
    n_medium : float
        Refractive index of the medium (tissue), typically 1.4
    n_probe : ndarray, shape (N,)
        Refractive index at the probe (interpolated over lambda_wl)

    Returns
    -------
    res : ndarray, shape (N,)
        Computed reflectance
    """
    x = np.asarray(x, dtype=float)
    c1, c2, d = x[0], x[1], x[2]
    c_blood, StO2 = x[3], x[4]
    cyt_oxi, cyt_red = x[6], x[7]

    lam = np.asarray(lambda_wl, dtype=float)
    n_pr = np.asarray(n_probe, dtype=float)

    # Fresnel reflection coefficient
    n = n_medium / n_pr
    Rf = 0.0636 * n + 0.668 + 0.710 / n - 1.44 / n**2
    kd = (1.0 + Rf) / (1.0 - Rf)

    # Reduced scattering (modified power law)
    uscat = c1 * (d * (lam / 600.0)**(-4) + (1.0 - d) * (lam / 600.0)**(-c2))
    uscat = np.maximum(uscat, 1e-10)  # avoid division by zero

    # Blood absorption
    ublood = c_blood * (StO2 * ua[0, :] + (1.0 - StO2) * ua[1, :])

    # Total absorption (mua_fat not included, as in the commented-out MATLAB code)
    uabs = ublood + cyt_oxi * ua[3, :] + cyt_red * ua[4, :]
    uabs = np.maximum(uabs, 1e-10)

    # Diffusion parameters
    a = uscat / (uscat + uabs)          # single scattering albedo
    z0 = 1.0 / uscat                    # virtual source depth
    D = 1.0 / (3.0 * uscat)            # diffusion coefficient
    zb = 2.0 * kd * D                  # extrapolated boundary thickness

    # Distances to image sources
    r1 = np.sqrt(z0**2 + rho**2)
    r2 = np.sqrt((z0 + 2.0 * zb)**2 + rho**2)

    # Effective attenuation coefficient
    ueff = np.sqrt(3.0 * uabs * uscat)

    # Two-source model (surface reflectance)
    res1 = z0 * (ueff + 1.0 / r1) * np.exp(-ueff * r1) / r1**2
    res2 = (z0 + 2.0 * zb) * (ueff + 1.0 / r2) * np.exp(-ueff * r2) / r2**2

    res = (a / (4.0 * np.pi)) * (res1 + res2)

    return res
