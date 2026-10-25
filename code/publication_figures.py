"""Rebuild manuscript/report figures from the verified validation checkpoint.

Run after bvcm.ipynb and check_bioskin_rerun.py. No model is refitted.
PDFs retain vector text/axes; dense data layers are rasterized at 600 dpi.
"""
import hashlib
import argparse
import json
import os
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'output/corrected_validation'
MANUSCRIPT = ROOT.parent / 'SkinSpectra_manuscript'
os.environ.setdefault('MPLCONFIGDIR', str(ROOT / '.runtime/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
from scipy.special import expit
from scipy.stats import chi2
from sklearn.mixture import GaussianMixture
from bioskin_inputs import load_corrected_bioskin
from check_bioskin_rerun import evaluate

BLUE, ORANGE, GREEN = '#0072B2', '#D55E00', '#009E73'
PALETTE = [BLUE, ORANGE, GREEN, '#CC79A7', '#8064A2', '#555555']
plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 8.5, 'axes.titlesize': 9,
    'axes.labelsize': 8.5, 'legend.fontsize': 7, 'xtick.labelsize': 8,
    'ytick.labelsize': 8, 'axes.spines.top': False, 'axes.spines.right': False,
    'axes.linewidth': .6, 'lines.linewidth': 1.1, 'xtick.major.width': .6,
    'ytick.major.width': .6, 'xtick.major.size': 3, 'ytick.major.size': 3,
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'savefig.dpi': 600,
    'legend.frameon': False, 'figure.facecolor': 'white', 'axes.facecolor': 'white',
})
# Match the current NeurIPS manuscript text width to preserve printed font sizes.
WIDTH = 5.5
REPORT_WIDTH = 180 / 25.4
manifest = []


def write_manifest(coverage):
    path = OUT / 'figure_manifest.json'
    # Preserve entries for informal report figures when only the paper is updated.
    previous = json.loads(path.read_text())['exports'] if path.exists() else []
    changed = {entry['file'] for entry in manifest}
    exports = [e for e in previous if e['file'] not in changed] + manifest
    path.write_text(json.dumps({'dpi':600, 'font':'DejaVu Sans', 'exports':exports,
        'pca_coverage_percent':coverage,
        'reflectance_colormap':{'colors':['#fbf7ec','#a15c33','#42230c'],
                                'mi_limits':[20,150], 'clip':True}}, indent=2))


def spectrum_colorbar(fig, axes, norm, cmap):
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axes,
                      orientation='horizontal', fraction=.05, aspect=45,
                      shrink=.85, pad=.06, label='Melanin index')
    bar.set_ticks([20, 50, 75, 100, 125, 150])


def save(fig, name, manuscript=True):
    fig.canvas.draw()
    # Catch cropped annotations before writing publication files.
    frame = fig.bbox
    artists = list(fig.texts)
    for ax in fig.axes:
        artists += [ax.title, ax._left_title, ax._right_title,
                    ax.xaxis.label, ax.yaxis.label, *ax.texts]
        if ax.get_legend():
            artists += [ax.get_legend().get_title(), *ax.get_legend().get_texts()]
    for artist in artists:
        if artist.get_visible() and artist.get_text():
            box = artist.get_window_extent(fig.canvas.get_renderer())
            if box.width and box.height:
                assert (box.x0 >= frame.x0-1 and box.y0 >= frame.y0-1
                        and box.x1 <= frame.x1+1 and box.y1 <= frame.y1+1), (name, artist.get_text())
    folder = ROOT / 'figures' if manuscript else OUT / 'figures'
    folder.mkdir(parents=True, exist_ok=True)
    for ext in ['pdf', 'png']:
        path = folder / f'{name}.{ext}'
        fig.savefig(path, dpi=600, facecolor='white')
        # The manuscript includes PDFs only; PNGs stay in figures/ for the website.
        if manuscript and ext == 'pdf' and (MANUSCRIPT / 'figures').is_dir():
            shutil.copy2(path, MANUSCRIPT / 'figures' / path.name)
        manifest.append({'file': str(path.relative_to(ROOT)),
                         'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                         'width_inches': float(fig.get_figwidth()),
                         'height_inches': float(fig.get_figheight())})
    plt.close(fig)


def panel(ax, letter, title):
    ax.set_title(f'({letter})  {title}', loc='left', pad=8)


def ellipse(ax, mean, cov, **kwargs):
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    angle = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    wh = 2 * np.sqrt(chi2.ppf(.95, 2) * vals)
    ax.add_patch(Ellipse(mean, *wh, angle=angle, fill=False, **kwargs))


def spectra_layer(ax, wavelengths, spectra, mi, norm, cmap):
    segments = np.stack([np.broadcast_to(wavelengths, spectra.shape), spectra], axis=-1)
    lines = LineCollection(segments, array=mi, norm=norm, cmap=cmap,
                           linewidths=.35, alpha=.70, rasterized=True)
    ax.add_collection(lines)
    ax.set(xlim=(400, 700), xlabel='Wavelength (nm)')
    ax.set_xticks([400, 500, 600, 700])
    return lines


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manuscript-only', action='store_true')
    args = parser.parse_args()
    if not (MANUSCRIPT / 'figures').is_dir():
        print(f'Warning: {MANUSCRIPT / "figures"} not found; manuscript figures will not be copied.')
    summary =json.loads((OUT / 'summary.json').read_text())
    ref = dict(np.load(OUT / 'reference_state.npz'))
    checkpoint = dict(np.load(ROOT / 'data/bvcm_gibbs_validation_checkpoint.npz'))
    # Section 3 figures use the training-split fit; post_betas_mean is the full-data validation reference.
    wl, beta = ref['wavelengths'], ref['evaluation_post_betas_mean']
    workbook = ROOT / 'data/ISSA_17_Jan_2025_Yan_Lu.xlsx'
    assert hashlib.sha256(workbook.read_bytes()).hexdigest() == summary['issa_workbook_sha256']
    data = pd.read_excel(workbook, sheet_name='ISSA')
    sites = data.iloc[11:, 6].to_numpy(dtype=int)
    face = np.isin(sites, [2, 3, 4, 5, 6, 9])
    spectra = data.iloc[11:, 17:48].to_numpy(dtype=float)[face]
    sites = sites[face]
    derivative = savgol_filter(spectra, 5, 3, axis=1, deriv=1) / np.gradient(wl)
    keep = ~np.any(abs(derivative - derivative.mean(0)) > 5 * derivative.std(0), axis=1)
    spectra, sites = spectra[keep], sites[keep]
    assert len(spectra) == 7923
    mi = 100 * np.log10(100 / spectra[:, 28])
    ei = 100 * np.log10(spectra[:, 28] / spectra[:, 16])
    corrected = load_corrected_bioskin(ROOT / 'bioskin_data', wl)
    assert corrected['source_files'] == summary['source_files']
    new = corrected['spectra_percent']
    old_full = np.loadtxt(ROOT / 'new_sample_dataset_package/new_sample_spectral_reflectance_smoothed.csv', delimiter=',')
    old = 100 * old_full[:, np.searchsorted(380 + 2 * np.arange(310), wl)]
    a = evaluate(old, ref, summary['stage1_threshold'])
    b = evaluate(new, ref, summary['stage1_threshold'])
    # Counts must match the current (full-data) Stage 2 reference exported by bvcm.ipynb.
    assert b['stage1'].sum() == summary['stage1_pass'] and b['final'].sum() == summary['passed_both']
    # Preserve the original ISSA/bvcm notebook palette and MI normalization.
    cmap = LinearSegmentedColormap.from_list('fitzpatrick', ['#fbf7ec', '#a15c33', '#42230c'])
    norm = Normalize(20, 150, clip=True)

    fig, axes = plt.subplots(2, 3, figsize=(WIDTH, 3.6), sharex=True, sharey=True, layout='constrained')
    names = ['Cheek', 'Cheekbone', 'Chin', 'Earlobe', 'Forehead', 'Nose tip']
    for ax, site, name, letter, color in zip(axes.flat, [2, 3, 4, 5, 6, 9], names, 'abcdef', PALETTE):
        mask = sites == site
        ax.scatter(mi[mask], ei[mask], s=3, alpha=.45, color=color, edgecolors='none', rasterized=True)
        panel(ax, letter, f'{name}\nn = {mask.sum():,}')
        ax.set(xlim=(0, 125), ylim=(0, 65))
    for ax in axes[-1]: ax.set_xlabel('Melanin index')
    for ax in axes[:, 0]: ax.set_ylabel('Erythema index')
    save(fig, 'melanin_vs_erythema_index_on_face_sites')

    fig, ax = plt.subplots(figsize=(WIDTH, 2.7), layout='constrained')
    ax.scatter(mi, ei, s=2, color='#888888', alpha=.25, edgecolors='none', rasterized=True)
    for i, (mean, cov, color) in enumerate(zip(ref['gmm_means'], ref['gmm_covariances'], PALETTE)):
        ellipse(ax, mean, cov, edgecolor=color, lw=1.2, label=f'Component {i+1}')
        ax.scatter(*mean, color=color, marker='+', s=36, linewidths=1.2, zorder=4)
    ax.set(xlabel='Melanin index', ylabel='Erythema index', xlim=(0, 125), ylim=(0, 65))
    ax.legend(loc='center left', bbox_to_anchor=(1, .5), title='95% component ellipses', title_fontsize=8)
    save(fig, 'gmm_mi_ei_confidence_regions')

    post = checkpoint['beta_curves']
    assert np.allclose(post.mean(0), beta)
    low, high = np.percentile(post, [2.5, 97.5], axis=0)
    fig, axes = plt.subplots(1, 3, figsize=(WIDTH, 2.25), layout='constrained')
    for k, (ax, title, color) in enumerate(zip(axes, [r'Intercept $\beta_0$', r'Melanin $\beta_{\rm MI}$', r'Erythema $\beta_{\rm EI}$'], PALETTE)):
        ax.fill_between(wl, low[k], high[k], color=color, alpha=.22)
        ax.plot(wl, beta[k], color=color)
        ax.set(xlabel='Wavelength (nm)', xlim=(400, 700))
        ax.set_xticks([400, 500, 600, 700])
        panel(ax, 'abc'[k], title)
    axes[0].set_ylabel('Coefficient (logit scale)')
    save(fig, 'beta_curves_estimation')

    fig, axes = plt.subplots(1, 3, figsize=(WIDTH, 2.6), sharey=True, layout='constrained')
    predicted_max = 0
    for ax, q, letter in zip(axes, [5, 50, 95], 'abc'):
        e = np.percentile(ei, q)
        for m in np.percentile(mi, [1, 5, 25, 50, 75, 95, 99]):
            curve = 100 * expit(np.array([1, m, e]) @ beta)
            predicted_mi = 100 * np.log10(100 / curve[28])
            ax.plot(wl, curve, color=cmap(norm(predicted_mi)), lw=1.3)
            predicted_max = max(predicted_max, curve.max())
        panel(ax, letter, f'EI: {q}th percentile')
        ax.set(xlabel='Wavelength (nm)', xlim=(400, 700))
        ax.set_xticks([400, 500, 600, 700])
    axes[0].set_ylabel('Reflectance (%)')
    axes[0].set_ylim(0, np.ceil(predicted_max/5)*5)
    spectrum_colorbar(fig, axes, norm, cmap)
    save(fig, 'bvcm_reflectance_predictions_for_different_mi_percentiles_given_ei')

    test = ref['test_indices']
    predicted = 100 * expit(np.column_stack([np.ones(len(test)), mi[test], ei[test]]) @ beta)
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH, 2.75), sharex=True, sharey=True, layout='constrained')
    ymax = np.ceil(max(predicted.max(), spectra[test].max()) / 10) * 10
    for ax, values, title, letter in zip(axes, [predicted, spectra[test]], ['BVCM prediction\nn = 1,585', 'ISSA observations\nn = 1,585'], 'ab'):
        spectra_layer(ax, wl, values, mi[test], norm, cmap)
        ax.set_ylim(0, ymax)
        panel(ax, letter, title)
    axes[0].set_ylabel('Reflectance (%)')
    spectrum_colorbar(fig, axes, norm, cmap)
    save(fig, 'predicted_vs_issa_test_spectra')

    physical = corrected['physical_parameters']
    columns = ['Melanin', 'Haemoglobin', 'Thickness_cm', 'MelaninTypeRatio', 'BloodOxygenation']
    titles = ['Melanin\nvolume fraction', 'Blood\nvolume fraction', 'Epidermal\nthickness', 'Melanin type\nratio', 'Blood\noxygenation']
    labels = ['Volume fraction', 'Volume fraction', 'Thickness (µm)', 'Fraction', 'Fraction']
    fig = plt.figure(figsize=(WIDTH, 3.65), layout='constrained')
    grid = fig.add_gridspec(2, 6)
    slots = [grid[0,0:2],grid[0,2:4],grid[0,4:6],grid[1,1:3],grid[1,3:5]]
    axes = [fig.add_subplot(slot) for slot in slots]
    for k, (ax, col, title, xlabel) in enumerate(zip(axes, columns, titles, labels)):
        values = physical[col].to_numpy() * (10000 if col == 'Thickness_cm' else 1)
        ax.hist(values, bins=20, color=BLUE, alpha=.85, edgecolor='white', linewidth=.35)
        panel(ax, 'abcde'[k], title)
        ax.set(xlabel=xlabel, ylim=(0,260), yticks=[0,100,200])
        if k in [0,3]: ax.set_ylabel('Frequency')
        else: ax.tick_params(labelleft=False)
    save(fig, 'bioskin_properties_distribution')

    fig, axes = plt.subplots(1, 3, figsize=(WIDTH, 2.65), sharex=True, sharey=True, layout='constrained')
    masks = [np.ones(len(new), dtype=bool), b['final'], ~b['final']]
    ymax = np.ceil(new.max()/10)*10
    for ax, mask, title, letter in zip(axes, masks, ['All candidates', 'Retained', 'Rejected'], 'abc'):
        spectra_layer(ax, wl, new[mask], b['mi'][mask], norm, cmap)
        ax.set_ylim(0, ymax)
        panel(ax, letter, f'{title}\nn = {mask.sum():,}')
    axes[0].set_ylabel('Reflectance (%)')
    spectrum_colorbar(fig, axes, norm, cmap)
    save(fig, 'bioskin_spectra_validation_overview')

    gmm = GaussianMixture(n_components=5, covariance_type='full')
    gmm.weights_, gmm.means_, gmm.covariances_ = ref['gmm_weights'], ref['gmm_means'], ref['gmm_covariances']
    np.random.seed(42)
    sampled, _ = gmm.sample(3000)
    synthetic = 100 * expit(np.column_stack([np.ones(3000), sampled]) @ beta)
    scores = [(s - ref['pca_mean']) @ ref['pca_components'].T for s in [spectra, synthetic, new]]
    coverage = []
    for s in scores:
        centered = s - ref['pca_score_mean']
        coverage.append(100 * np.mean(np.einsum('ij,jk,ik->i', centered, np.linalg.inv(ref['pca_score_cov']), centered) <= chi2.ppf(.95, 2)))
    assert np.allclose(coverage, [93.39896503849552, 95.4, 22.775])
    variance = np.var(scores[0], axis=0, ddof=1) / np.var(spectra, axis=0, ddof=1).sum() * 100
    all_scores = np.concatenate(scores)
    lo, hi = all_scores.min(0), all_scores.max(0)
    margin = .06 * (hi-lo)
    fig, axes = plt.subplots(1, 3, figsize=(WIDTH, 2.45), sharex=True, sharey=True, layout='constrained')
    for i, (ax, s, title, color) in enumerate(zip(axes, scores, ['ISSA\nn = 7,923', 'BVCM\nn = 3,000', 'BioSkin\nn = 4,000'], [BLUE, GREEN, ORANGE])):
        if i:
            ax.scatter(*scores[0].T, s=2, c='.8', alpha=.25, edgecolors='none', rasterized=True)
        ax.scatter(*s.T, s=3, color=color, alpha=.35, edgecolors='none', rasterized=True)
        ellipse(ax, ref['pca_score_mean'], ref['pca_score_cov'], edgecolor='black', lw=.9, ls='--', zorder=4)
        panel(ax, 'abc'[i], title)
        ax.set(xlim=(lo[0]-margin[0], hi[0]+margin[0]), ylim=(lo[1]-margin[1], hi[1]+margin[1]), xlabel=f'PC1 ({variance[0]:.1f}%)')
        ax.text(.04, .05, f'Within ellipse\n{coverage[i]:.1f}%', transform=ax.transAxes, fontsize=7.5, bbox=dict(facecolor='white', edgecolor='none', alpha=.85, pad=2))
        ax.set_xticks([-150, 0, 150])
    axes[0].set_ylabel(f'PC2 ({variance[1]:.1f}%)')
    save(fig, 'pca_issa_bvcm_bioskin_comparison')

    if args.manuscript_only:
        write_manifest(coverage)
        print('Exported 8 manuscript figures; original MI palette restored; counts and PCA coverage verified.')
        return

    fig, axes = plt.subplots(1, 2, figsize=(REPORT_WIDTH, 3.15), gridspec_kw={'width_ratios':[1, 1.2]}, layout='constrained')
    counts = np.array([[734, 3121, 145], [733, 2833, 434]])
    bottom = np.zeros(2)
    for j, (label, color) in enumerate(zip(['Retained', 'Rejected: Stage 1', 'Rejected: Stage 2'], [GREEN, '#B6BDC3', ORANGE])):
        axes[0].bar([0, 1], 100*counts[:, j]/4000, bottom=bottom, color=color, width=.57, label=label)
        for i in range(2):
            midpoint = bottom[i]+100*counts[i, j]/8000
            if counts[i,j] < 200:
                axes[0].annotate(f'{counts[i,j]:,}', xy=(i+.285,midpoint),
                                 xytext=(i+.49,midpoint), ha='center', va='center',
                                 fontsize=7, arrowprops={'arrowstyle':'-', 'lw':.5})
            else:
                axes[0].text(i, midpoint, f'{counts[i,j]:,}', ha='center', va='center', fontsize=8)
        bottom += 100*counts[:, j]/4000
    axes[0].set(xticks=[0,1], xticklabels=['Previous', 'Corrected'], ylabel='Percentage of 4,000 spectra', ylim=(0,100))
    axes[0].legend(loc='upper center', bbox_to_anchor=(.5,-.14), ncol=1)
    panel(axes[0], 'a', 'Overall outcomes')
    transition = pd.crosstab(a['outcome'], b['outcome']).reindex(index=['passed_both','rejected_stage_1','rejected_stage_2'], columns=['passed_both','rejected_stage_1','rejected_stage_2']).to_numpy()
    axes[1].imshow(transition, cmap='Blues', norm=matplotlib.colors.PowerNorm(.45), aspect='equal')
    for i in range(3):
        for j in range(3):
            axes[1].text(j,i,f'{transition[i,j]:,}',ha='center',va='center',color='white' if transition[i,j]>1200 else 'black',fontsize=10)
    labels = ['Retained','Rejected\nStage 1','Rejected\nStage 2']
    axes[1].set(xticks=range(3),yticks=range(3),xticklabels=labels,yticklabels=labels,xlabel='Corrected outcome',ylabel='Previous outcome')
    axes[1].tick_params(length=0)
    panel(axes[1], 'b', 'Paired transitions (counts)')
    save(fig, 'rerun_outcomes_and_transitions', manuscript=False)

    curves = pd.read_csv(OUT / 'changed_acceptance_example.csv')
    diag = pd.read_csv(OUT / 'lost_stage2_wavelength_diagnostics.csv')
    fig, axes = plt.subplots(2, 2, figsize=(REPORT_WIDTH, 5.0), layout='constrained')
    for prefix, color, label in [('old',BLUE,'Previous'),('new',ORANGE,'Corrected')]:
        axes[0,0].plot(wl,curves[f'{prefix}_reflectance_percent'],color=color,label=f'{label}: spectrum')
        axes[0,0].plot(wl,curves[f'{prefix}_bvcm_prediction_percent'],color=color,ls='--',label=f'{label}: prediction')
        axes[0,1].plot(wl,curves[f'{prefix}_logit_residual'],color=color,label=label)
        axes[1,0].plot(wl,diag[f'{prefix}_conditional_z_rms'],color=color,label=label)
    axes[0,0].legend(loc='upper left',fontsize=6.7)
    axes[0,0].set_ylabel('Reflectance (%)')
    axes[0,1].set_ylabel('Residual (logit scale)')
    axes[0,1].axhline(0,color='.65',lw=.6)
    axes[0,1].legend()
    axes[1,0].set_ylabel('RMS conditional standardized residual')
    axes[1,0].legend()
    panel(axes[0,0],'a','Sample 2408: spectrum and prediction')
    panel(axes[0,1],'b','Sample 2408: conditional residual')
    panel(axes[1,0],'c','184 spectra newly rejected at Stage 2')
    for ax in [axes[0,0],axes[0,1],axes[1,0]]:
        ax.set(xlim=(400,700),xlabel='Wavelength (nm)')
    mask = a['final'] & b['stage1']
    for group,color,label in [(mask & b['final'],'.55','Retained in both (539)'),(mask & ~b['final'],ORANGE,'New Stage 2 rejection (184)')]:
        axes[1,1].scatter(a['spectral_distance'][group],b['spectral_distance'][group],s=7,color=color,alpha=.65,edgecolors='none',label=label,rasterized=True)
    cutoff = summary['stage2_threshold']
    axes[1,1].axhline(cutoff,color='black',ls='--',lw=.8)
    axes[1,1].axvline(cutoff,color='black',ls='--',lw=.8)
    axes[1,1].plot([0,cutoff],[0,cutoff],ls=':',color='.5',lw=.7)
    axes[1,1].set(xlim=(0,cutoff*1.05),ylim=(0,max(16,b['spectral_distance'][mask].max()*1.05)),xlabel='Previous Stage 2 distance',ylabel='Corrected Stage 2 distance')
    axes[1,1].legend(loc='upper left',fontsize=6.7)
    panel(axes[1,1],'d','Paired spectral distances')
    save(fig,'rerun_changed_acceptance',manuscript=False)
    write_manifest(coverage)
    print(f'Exported {len(manifest)//2} figures as PDF and 600-dpi PNG; verified counts, PCA coverage and input hashes.')


if __name__ == '__main__':
    main()
