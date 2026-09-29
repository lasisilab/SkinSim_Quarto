"""Export the coauthor Markdown report to self-contained HTML.

Word export requires an explicit --docx flag. Does not invoke TeX or create PDFs.
"""
import argparse
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pandoc', type=Path)
    parser.add_argument('--docx', action='store_true', help='Also export Word, only when explicitly requested.')
    args = parser.parse_args()
    pandoc = args.pandoc or shutil.which('pandoc')
    if not pandoc:
        quarto = shutil.which('quarto')
        candidates = [Path(quarto).parent / 'tools/pandoc.exe'] if quarto else []
        candidates += [Path('C:/Program Files/Quarto/bin/tools/pandoc.exe')]
        pandoc = next((p for p in candidates if p.is_file()), None)
    if not pandoc:
        raise SystemExit('Pandoc is not available; use --pandoc to supply its path.')
    out = ROOT / 'output/corrected_validation'
    common = [str(pandoc), 'rerun_report.md', '--from=markdown', '--standalone']
    formats = [
        ['--to=html5', '--embed-resources', '--css=report_style.css', '--output=rerun_report.html'],
    ]
    if args.docx:
        formats.append(['--to=docx', '--output=rerun_report.docx'])
    for extra in formats:
        subprocess.run(common + extra, cwd=out, check=True)
    print('Exported rerun_report.html (embedded figures)' + (' and rerun_report.docx' if args.docx else '') + '. No TeX/PDF conversion.')


if __name__ == '__main__':
    main()
