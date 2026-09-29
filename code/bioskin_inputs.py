"""Read and check Carlos's corrected BioSkin CSVs before validation."""

from hashlib import sha256
from pathlib import Path

import numpy as np
import pandas as pd


def load_corrected_bioskin(directory, wavelengths):
    """Return aligned reflectance (%) and both parameter representations.

    Network columns are melanin, haemoglobin, normalized thickness, melanin
    type ratio, blood oxygenation, and occlusion. The physical-units CSV has
    its own header order, with oxygenation before melanin type ratio.
    """
    directory = Path(directory)
    paths = {
        "spectra": directory / "corrected_spectral_reflectance_smoothed.csv",
        "network_parameters": directory / "corrected_skin_props_network_ready.csv",
        "physical_parameters": directory / "corrected_skin_props_raw_physical_units.csv",
    }
    spectra = np.loadtxt(paths["spectra"], delimiter=",")
    network = np.loadtxt(paths["network_parameters"], delimiter=",")
    physical = pd.read_csv(paths["physical_parameters"])
    physical_columns = [
        "Melanin", "Haemoglobin", "Thickness_cm",
        "BloodOxygenation", "MelaninTypeRatio",
    ]
    if list(physical.columns) != physical_columns:
        raise ValueError("Unexpected physical-parameter headers in corrected CSV.")
    if spectra.shape != (4000, 310) or network.shape != (4000, 6):
        raise ValueError("Expected 4,000 corrected spectra (310 wavelengths) and six network inputs.")
    if physical.shape != (4000, 5):
        raise ValueError("Physical parameters must have one row per corrected spectrum.")
    if not all(np.isfinite(a).all() for a in (spectra, network, physical.to_numpy())):
        raise ValueError("BioSkin inputs contain non-finite values.")
    if not ((spectra > 0) & (spectra < 1)).all():
        raise ValueError("Corrected reflectance must be in (0, 1) before percentage conversion.")

    # CSV rounding differs slightly across the supplied parameter files.
    # These checks catch swapped slots and incorrect thickness scaling.
    expected_physical = np.column_stack([
        network[:, 0] ** 3,
        network[:, 1] ** 4,
        0.001 + 0.034 * network[:, 2],
        network[:, 4],
        network[:, 3],
    ])
    if not np.allclose(expected_physical, physical.to_numpy(), rtol=0, atol=1e-7):
        raise ValueError("Network and physical parameters disagree; check slot order and scaling.")
    if not np.allclose(network[:, 5], 1.0, rtol=0, atol=1e-8):
        raise ValueError("Expected occlusion=1 in the supplied network-ready inputs.")
    bounds = {
        "Thickness_cm": (0.005, 0.025),
        "MelaninTypeRatio": (0.65, 0.80),
        "BloodOxygenation": (0.90, 1.00),
    }
    for column, (low, high) in bounds.items():
        if not physical[column].between(low - 1e-7, high + 1e-7).all():
            raise ValueError(f"{column} is outside the corrected dataset's stated range.")

    source_wavelengths = 380.0 + 2 * np.arange(spectra.shape[1])
    wavelengths = np.asarray(wavelengths, dtype=float)
    if not np.array_equal(wavelengths, np.arange(400.0, 701.0, 10.0)):
        raise ValueError("This validation expects the 31 ISSA wavelengths from 400 to 700 nm.")
    columns = np.searchsorted(source_wavelengths, wavelengths)
    if not np.array_equal(source_wavelengths[columns], wavelengths):
        raise ValueError("BioSkin and ISSA wavelengths do not align exactly.")

    return {
        "spectra_percent": 100.0 * spectra[:, columns],
        "network_parameters": network,
        "physical_parameters": physical,
        "source_wavelengths": source_wavelengths,
        "source_files": {
            key: {"name": path.name, "sha256": sha256(path.read_bytes()).hexdigest()}
            for key, path in paths.items()
        },
    }
