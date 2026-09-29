# SkinSim

SkinSim is a Quarto-based research website for statistical modeling of human skin reflectance spectra.

The project integrates real measurements (ISSA) and simulated spectra to support:

- reproducible spectral analysis,
- biologically informed modeling (melanin and erythema effects),
- synthetic spectrum generation and validation.

## Project structure

- `analysis/`: Quarto analysis reports (`.qmd`)
- `notebooks/`: Jupyter notebooks (`.ipynb`)
- `code/`: reusable Python and R modules/scripts
- `data/`: source and derived datasets
- `output/`: generated figures, tables, and exports
- `docs/`: rendered static website output

## Environment setup

### R environment

```r
source("renv-setup.R")
```

### Python environment

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Local development

Run the live preview server:

```bash
quarto preview
```

## Render website

Published site: <https://lasisilab.github.io/SkinSim_Quarto/>

- [BioSkin validation](https://lasisilab.github.io/SkinSim_Quarto/analysis/bioskin_validation.html)
- [BVCM notebook](https://lasisilab.github.io/SkinSim_Quarto/notebooks/bvcm.html)
- [Rerun report](https://lasisilab.github.io/SkinSim_Quarto/analysis/rerun_report.html)

To publish already-executed notebook outputs without refitting the model,
render with `quarto render --no-execute --no-clean`. For partial updates,
pass the changed notebook or `.qmd` path to that command. Commit the updated
sources and `docs/` output; GitHub Pages serves the generated HTML.

Regenerate the standalone report with
`python code/export_rerun_report.py` (HTML only). Word output requires an
explicit `--docx` flag. The raw ISSA workbook, local environments, temporary
files, and model checkpoints in `data/` are excluded from Git.

Build the full site to `docs/`:

```bash
quarto render
```

## Recommended workflow

1. Add or update analysis content in `analysis/` or `notebooks/`.
2. Keep reusable model code in `code/`.
3. Store intermediate and final artifacts in `output/`.
4. Preview frequently with `quarto preview`.
5. Commit both source files and rendered updates when publishing.

## Notes

- The website uses Quarto navigation and supports both R and Python execution.
- Notebook pages are rendered as part of the same site output.
