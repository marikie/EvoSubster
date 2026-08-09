#!/bin/bash

set -eu

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/evosubster-one2one.XXXXXX")
trap 'rm -rf "$tmp_dir"' EXIT

stub_bin="$tmp_dir/bin"
mkdir -p "$stub_bin" "$tmp_dir/db"

cat > "$stub_bin/lastal" <<'EOF'
#!/bin/bash
printf '%s\n' "$@" > "$LASTAL_ARGS_LOG"
printf '##maf version=1\n\n'
EOF

cat > "$stub_bin/last-split" <<'EOF'
#!/bin/bash
cat
EOF

cat > "$stub_bin/maf-linked" <<'EOF'
#!/bin/bash
cat "$1"
EOF

cat > "$stub_bin/lastdb" <<'EOF'
#!/bin/bash
exit 0
EOF

cat > "$stub_bin/last-train" <<'EOF'
#!/bin/bash
exit 0
EOF

chmod +x "$stub_bin/lastal" "$stub_bin/last-split" "$stub_bin/maf-linked" \
    "$stub_bin/lastdb" "$stub_bin/last-train"

reference="$tmp_dir/reference.fna"
query="$tmp_dir/query.fna"
train="$tmp_dir/params.train"
output="$tmp_dir/output.maf"
printf '>reference\nACGT\n' > "$reference"
printf '>query\nACGT\n' > "$query"
printf '# parameters\n' > "$train"

PATH="$stub_bin:$PATH" LASTAL_ARGS_LOG="$tmp_dir/lastal.args" THREAD_NUM=4 \
    bash "$ROOT_DIR/src/align/one2one.sh" 20260809 "$reference" "$query" \
    "$tmp_dir/db" "$train" "$output" > "$tmp_dir/run.log" 2>&1

test -s "$output"
grep -Fxq -- '-i' "$tmp_dir/lastal.args"
grep -Fxq -- '1M' "$tmp_dir/lastal.args"
grep -Fxq -- '--split-b=8G' "$tmp_dir/lastal.args"
test "$(grep -Fxc -- '-P4' "$tmp_dir/lastal.args")" -eq 1

echo "ok: one2one bounds lastal query batches while preserving caller thread control"
