#!/usr/bin/env python3
"""Generate a complete plot set privately; publish only successful jobs."""
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

from run_manifest import atomic_json, sha256


def generate(args):
    if len(args) != 11:
        raise ValueError('Expected eight TSVs, R directory, and two figure directories')
    scripts = Path(args[8]).resolve()
    font_script = scripts / 'required_fonts.R'
    subprocess.run(['Rscript', '-e', 'source(commandArgs(TRUE)[1]); load_required_fonts()', str(font_script)], check=True)
    inventory = []
    with tempfile.TemporaryDirectory(prefix='.plot-staging-', dir=Path.cwd()) as staging:
        stage = Path(staging)
        for offset, slot in ((0, 'org2'), (1, 'org3')):
            output = Path(args[9 + offset])
            for noncoding in (False, True):
                index = offset + (2 if noncoding else 0)
                single = Path(args[index]).resolve()
                double = Path(args[index + 4]).resolve()
                if noncoding and not single.exists() and not double.exists():
                    continue
                if not single.is_file() or not double.is_file():
                    raise ValueError('Missing paired single/dinucleotide plot inputs')
                work = stage / slot / ('ncds' if noncoding else 'all')
                single_out = work / 'singlenuc'
                log_out = single_out / 'log-ratio'
                dinuc_out = work / 'dinuc'
                for path in (single_out, log_out, dinuc_out):
                    path.mkdir(parents=True, exist_ok=True)
                for script, source, dest in (('sbmut.R', single, single_out),
                                              ('logRatioPlot.R', single, log_out),
                                              ('dinucleotide-plot.R', double, dinuc_out)):
                    subprocess.run(['Rscript', str(scripts / script), str(source), str(dest)], check=True)
                suffix = 'ncds_' if noncoding else ''
                expected = [(single_out / 'ratio' / (single.stem + '_norm.pdf'), 'singlenuc/ratio', 'pdf_' + suffix + 'norm'),
                            (single_out / 'count' / (single.stem + '_sbst.pdf'), 'singlenuc/count', 'pdf_' + suffix + 'sbst'),
                            (single_out / 'count' / (single.stem + '_ori.pdf'), 'singlenuc/count', 'pdf_' + suffix + 'ori'),
                            (log_out / (single.stem + '_logRatio.pdf'), 'singlenuc/log-ratio', 'pdf_' + suffix + 'logratio'),
                            (dinuc_out / (double.name + '.pdf'), 'dinuc', 'pdf_dinuc_' + ('ncds_tsv' if noncoding else 'tsv'))]
                for source, subdirectory, kind in expected:
                    content = source.read_bytes()
                    if not content.startswith(b'%PDF-') or b'%%EOF' not in content[-1024:]:
                        raise ValueError(f'Incomplete PDF: {source}')
                    inventory.append((source, output / subdirectory / source.name, kind, slot))
        records = []
        for source, target, kind, slot in inventory:
            target.parent.mkdir(parents=True, exist_ok=True)
            source.replace(target)
            records.append(dict(path=target.as_posix(), kind=kind, slot=slot, sha256=sha256(target)))
        atomic_json(Path('metadata/plot_completion.json'), dict(status='complete', artifacts=records))


if __name__ == '__main__':
    generate(sys.argv[1:])
