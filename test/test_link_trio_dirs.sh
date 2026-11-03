#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
python3 - "$ROOT_DIR" <<'PY'
from pathlib import Path
import subprocess
import sys
import tempfile
root=Path(sys.argv[1])
sys.path.insert(0,str(root/'test'))
from test_collect_run_summary import create_run
with tempfile.TemporaryDirectory(dir=root/'tmp/codex',prefix='manifest-links-') as tmp:
    work=Path(tmp)
    first=create_run(work/'lineage'/'arbitrary')
    output=work/'links'
    command=['bash',str(root/'src/report/link_trio_dirs.sh'),'--output-dir',str(output),str(work/'lineage')]
    subprocess.run(command,check=True)
    link=output/'fixture-trio'
    assert link.is_symlink() and link.resolve()==first.resolve()
    subprocess.run(command,check=True)
    create_run(work/'other'/'different-folder')
    conflict=subprocess.run(command+[str(work/'other')],capture_output=True)
    assert conflict.returncode != 0 and link.resolve()==first.resolve()
    help_result=subprocess.run(['bash',str(root/'src/report/link_trio_dirs.sh'),'--help'],capture_output=True)
    assert help_result.returncode==0
print('All manifest trio linking tests passed.')
PY
