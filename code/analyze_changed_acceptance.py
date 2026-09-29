"""Diagnose spectra retained before the decoder correction but rejected after it."""
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from bioskin_inputs import load_corrected_bioskin
from check_bioskin_rerun import evaluate


def main():
    root = Path(__file__).resolve().parents[1]
    out = root / 'output/corrected_validation'
    summary = json.loads((out / 'summary.json').read_text(encoding='utf-8'))
    with np.load(out / 'reference_state.npz', allow_pickle=False) as saved:
        ref = dict(saved)
    loaded = load_corrected_bioskin(root / 'bioskin_data', ref['wavelengths'])
    new = loaded['spectra_percent']
    old_full = np.loadtxt(root / 'new_sample_dataset_package/new_sample_spectral_reflectance_smoothed.csv', delimiter=',')
    columns = np.searchsorted(380 + 2 * np.arange(old_full.shape[1]), ref['wavelengths'])
    old = 100 * old_full[:, columns]
    a = evaluate(old, ref, summary['stage1_threshold'])
    b = evaluate(new, ref, summary['stage1_threshold'])
    masks = {
        'lost_at_stage1': a['final'] & ~b['stage1'],
        'lost_at_stage2': a['final'] & b['stage1'] & ~b['final'],
        'retained_in_both': a['final'] & b['final'],
    }
    loss = a['final'] & ~b['final']
    assert loss.sum() == 195

    def describe(x):
        return {key: float(value) for key, value in zip(
            ['min', 'p25', 'median', 'p75', 'max'], np.percentile(x, [0, 25, 50, 75, 100]),
        )}

    report = {'stage1_threshold': summary['stage1_threshold'], 'stage2_threshold': summary['stage2_threshold'], 'groups': {}}
    for name, mask in masks.items():
        group = {
            'n': int(mask.sum()),
            'old_MI': describe(a['mi'][mask]), 'new_MI': describe(b['mi'][mask]),
            'old_EI': describe(a['ei'][mask]), 'new_EI': describe(b['ei'][mask]),
            'delta_MI': describe((b['mi'] - a['mi'])[mask]),
            'delta_EI': describe((b['ei'] - a['ei'])[mask]),
            'old_stage1_distance': describe(a['stage1_distance'][mask]),
            'new_stage1_distance': describe(b['stage1_distance'][mask]),
            'old_stage2_distance': describe(a['spectral_distance'][mask]),
            'new_MI_ge_100': int((b['mi'][mask] >= 100).sum()),
        }
        if name != 'lost_at_stage1':
            group['new_stage2_distance'] = describe(b['spectral_distance'][mask])
            group['new_stage2_above_1_1_times_cutoff'] = int((b['spectral_distance'][mask] > 1.1 * summary['stage2_threshold']).sum())
        report['groups'][name] = group

    table = loaded['physical_parameters'].copy()
    table.insert(0, 'sample_id', np.arange(1, len(old) + 1))
    table['new_outcome'] = b['outcome']
    for prefix, values in [('old', a), ('new', b)]:
        for key in ['mi', 'ei', 'stage1_distance', 'spectral_distance']:
            table[f'{prefix}_{key}'] = values[key]
    table.loc[loss].to_csv(out / 'previously_retained_now_rejected.csv', index=False)

    residuals = []
    predictions = []
    for spectra, result in [(old, a), (new, b)]:
        pred = np.column_stack([np.ones(len(old)), result['mi'], result['ei']]) @ ref['post_betas_mean']
        predictions.append(pred)
        p = spectra / 100
        residuals.append(np.log(p / (1 - p)) - pred)
    precision = np.linalg.inv(ref['total_cov'])
    # Conditional standardized residual at each wavelength given all others.
    # These are diagnostic coordinates, not additive contributions to D^2.
    conditional_z = [r @ precision / np.sqrt(np.diag(precision)) for r in residuals]
    mask = masks['lost_at_stage2']
    wavelength_table = pd.DataFrame({'wavelength_nm': ref['wavelengths']})
    for prefix, r, z, spectra in zip(['old', 'new'], residuals, conditional_z, [old, new]):
        wavelength_table[f'{prefix}_median_logit_residual'] = np.median(r[mask], axis=0)
        wavelength_table[f'{prefix}_conditional_z_rms'] = np.sqrt(np.mean(z[mask] ** 2, axis=0))
        wavelength_table[f'{prefix}_median_reflectance_percent'] = np.median(spectra[mask], axis=0)
    wavelength_table['median_reflectance_change_percentage_points'] = np.median((new - old)[mask], axis=0)
    wavelength_table.to_csv(out / 'lost_stage2_wavelength_diagnostics.csv', index=False)
    top = wavelength_table.nlargest(6, 'new_conditional_z_rms')
    report['largest_new_conditional_residuals'] = top.to_dict(orient='records')

    # Choose an actual sample nearest the group's median corrected distance.
    indices = np.flatnonzero(mask)
    example = indices[np.argmin(np.abs(b['spectral_distance'][indices] - np.median(b['spectral_distance'][indices])))]
    report['example_sample_id'] = int(example + 1)
    report['example'] = table.iloc[example].to_dict()
    sample_curves = pd.DataFrame({'wavelength_nm': ref['wavelengths']})
    for prefix, spectra, pred, residual in zip(['old', 'new'], [old, new], predictions, residuals):
        sample_curves[f'{prefix}_reflectance_percent'] = spectra[example]
        sample_curves[f'{prefix}_bvcm_prediction_percent'] = 100 / (1 + np.exp(-pred[example]))
        sample_curves[f'{prefix}_logit_residual'] = residual[example]
    sample_curves.to_csv(out / 'changed_acceptance_example.csv', index=False)

    os.environ.setdefault('MPLCONFIGDIR', str(root / '.runtime/matplotlib'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), dpi=180)
    colors = ['#0072B2', '#D55E00']
    labels = ['Previous', 'Corrected']
    wl = ref['wavelengths']
    for prefix, color, label in zip(['old', 'new'], colors, labels):
        axes[0, 0].plot(wl, sample_curves[f'{prefix}_reflectance_percent'], color=color, label=f'{label} spectrum')
        axes[0, 0].plot(wl, sample_curves[f'{prefix}_bvcm_prediction_percent'], color=color, linestyle='--', label=f'{label} MI/EI prediction')
        axes[0, 1].plot(wl, sample_curves[f'{prefix}_logit_residual'], color=color, marker='.', label=label)
        axes[1, 0].plot(wl, wavelength_table[f'{prefix}_conditional_z_rms'], color=color, marker='.', label=label)
    axes[0, 0].set(title=f'Sample {example + 1}: spectrum and conditional prediction', ylabel='Reflectance (%)')
    axes[0, 0].legend(fontsize=8, loc='upper left')
    axes[0, 1].set(title=f'Sample {example + 1}: residual changes', ylabel='Observed minus predicted logit reflectance')
    axes[0, 1].axhline(0, color='0.6', linewidth=0.8)
    axes[0, 1].legend()
    axes[1, 0].set(title='Local shape diagnostics across the 184 Stage 2 losses', ylabel='RMS conditional standardized residual')
    axes[1, 0].legend()
    for ax in [axes[0, 0], axes[0, 1], axes[1, 0]]:
        ax.set_xlabel('Wavelength (nm)')
        ax.grid(alpha=0.2)
    ax = axes[1, 1]
    for name, color, label in [('retained_in_both', '#777777', 'Still retained (539)'), ('lost_at_stage2', '#D55E00', 'New Stage 2 rejection (184)')]:
        group_mask = masks[name]
        ax.scatter(a['spectral_distance'][group_mask], b['spectral_distance'][group_mask], s=12, color=color, alpha=0.6, label=label)
    cutoff = summary['stage2_threshold']
    ax.axhline(cutoff, color='black', linestyle='--', linewidth=1, label=f'Cutoff = {cutoff:.2f}')
    ax.axvline(cutoff, color='black', linestyle='--', linewidth=1)
    ax.plot([2, 15], [2, 15], color='0.65', linewidth=0.8, linestyle=':')
    ax.set(xlim=(2, 10), ylim=(2, 15.5), xlabel='Previous Stage 2 distance', ylabel='Corrected Stage 2 distance', title='Previously retained spectra that still pass Stage 1')
    ax.legend(fontsize=8, loc='upper left')
    ax.grid(alpha=0.2)
    fig.suptitle('Why previously retained spectra can fail after the decoder correction', fontsize=15)
    fig.text(0.5, 0.01, 'Same fitted BVCM, covariance and thresholds. MI/EI and predictions are recalculated for each version.\nConditional residual diagnostics are not additive contributions to the total distance.', ha='center', fontsize=9)
    fig.tight_layout(rect=(0, 0.055, 1, 0.96))
    fig.savefig(out / 'changed_acceptance_diagnostics.png', bbox_inches='tight')
    plt.close(fig)
    (out / 'changed_acceptance_diagnostics.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
