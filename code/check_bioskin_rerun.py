"""Verify exported validation and compare previous/corrected spectra on one fit.

Run after notebooks/bvcm.ipynb has completed. This uses its saved reference
state; it neither fits a new model nor changes either validation criterion.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular

from bioskin_inputs import load_corrected_bioskin


def evaluate(spectra, reference, stage1_threshold):
    """Independently recompute the two stages using Cholesky solves."""
    mi = 100 * np.log10(100 / spectra[:, 28])
    ei = 100 * np.log10(spectra[:, 28] / spectra[:, 16])
    covariates = np.column_stack((mi, ei))
    distances = []
    for mean, covariance in zip(reference['gmm_means'], reference['gmm_covariances']):
        standardized = solve_triangular(
            np.linalg.cholesky(covariance), (covariates - mean).T, lower=True,
        )
        distances.append(np.linalg.norm(standardized, axis=0))
    minimum_distance = np.min(distances, axis=0)
    stage1 = minimum_distance <= stage1_threshold
    predictions = np.column_stack((np.ones(len(spectra)), covariates)) @ reference['post_betas_mean']
    proportions = spectra / 100
    residuals = np.log(proportions / (1 - proportions)) - predictions
    standardized = solve_triangular(
        np.linalg.cholesky(reference['total_cov']), residuals[stage1].T, lower=True,
    )
    spectral_distance = np.full(len(spectra), np.nan)
    spectral_distance[stage1] = np.linalg.norm(standardized, axis=0)
    final = stage1 & (spectral_distance <= float(reference['mahal_threshold']))
    outcomes = np.where(~stage1, 'rejected_stage_1', np.where(final, 'passed_both', 'rejected_stage_2'))
    return {
        'mi': mi, 'ei': ei, 'stage1_distance': minimum_distance,
        'stage1': stage1, 'spectral_distance': spectral_distance,
        'final': final, 'outcome': outcomes,
    }


def main():
    root = Path(__file__).resolve().parents[1]
    out = root / 'output/corrected_validation'
    summary = json.loads((out / 'summary.json').read_text(encoding='utf-8'))
    exported = pd.read_csv(out / 'per_spectrum_results.csv')
    with np.load(out / 'reference_state.npz', allow_pickle=False) as saved:
        reference = dict(saved)
    corrected = load_corrected_bioskin(root / 'bioskin_data', reference['wavelengths'])
    assert corrected['source_files'] == summary['source_files']
    previous_dir = root / 'new_sample_dataset_package'
    previous = np.loadtxt(previous_dir / 'new_sample_spectral_reflectance_smoothed.csv', delimiter=',')
    previous_physical = pd.read_csv(previous_dir / 'new_sample_skin_props_raw_physical_units.csv')
    assert np.allclose(previous_physical, corrected['physical_parameters'], rtol=0, atol=1e-7)
    columns = np.searchsorted(380 + 2 * np.arange(previous.shape[1]), reference['wavelengths'])
    inputs = {'Previous (swapped slots)': 100 * previous[:, columns], 'Corrected': corrected['spectra_percent']}
    results = {name: evaluate(spectra, reference, summary['stage1_threshold']) for name, spectra in inputs.items()}
    current = results['Corrected']
    assert np.array_equal(exported.sample_id, np.arange(1, 4001))
    assert np.array_equal(exported.stage1_pass, current['stage1'])
    assert np.array_equal(exported.final_pass, current['final'])
    assert np.array_equal(exported.outcome, current['outcome'])
    assert np.array_equal(exported.stage2_evaluated, current['stage1'])
    assert exported.loc[~current['stage1'], 'stage2_pass'].isna().all()
    assert np.allclose(exported.MI, current['mi'], atol=1e-10)
    assert np.allclose(exported.EI, current['ei'], atol=1e-10)
    assert np.allclose(exported.stage1_min_mahalanobis, current['stage1_distance'], atol=1e-10)
    assert np.allclose(exported.stage2_mahalanobis, current['spectral_distance'], atol=1e-9, equal_nan=True)
    assert current['stage1'].sum() == summary['stage1_pass']
    assert current['final'].sum() == summary['passed_both']

    rows = []
    for name, result in results.items():
        rows.append({
            'dataset': name, 'n': len(result['mi']),
            'passed_stage1': int(result['stage1'].sum()),
            'rejected_stage1': int((~result['stage1']).sum()),
            'rejected_stage2': int((result['stage1'] & ~result['final']).sum()),
            'passed_both': int(result['final'].sum()),
            'passed_both_percent': 100 * result['final'].mean(),
            'retained_MI_min': float(result['mi'][result['final']].min()),
            'retained_MI_max': float(result['mi'][result['final']].max()),
            'retained_EI_min': float(result['ei'][result['final']].min()),
            'retained_EI_max': float(result['ei'][result['final']].max()),
        })
    comparison = pd.DataFrame(rows)
    comparison.to_csv(out / 'previous_vs_corrected.csv', index=False)
    previous_result = results['Previous (swapped slots)']
    paired = pd.DataFrame({
        'sample_id': np.arange(1, 4001),
        'previous_outcome': previous_result['outcome'],
        'corrected_outcome': current['outcome'],
        'previous_MI': previous_result['mi'], 'corrected_MI': current['mi'],
        'previous_EI': previous_result['ei'], 'corrected_EI': current['ei'],
    })
    paired.to_csv(out / 'paired_outcomes.csv', index=False)
    pd.crosstab(paired.previous_outcome, paired.corrected_outcome).to_csv(out / 'outcome_transitions.csv')
    print(comparison.to_string(index=False))
    print('Verified all 4,000 exported decisions and distances using independent Cholesky calculations.')


if __name__ == '__main__':
    main()
