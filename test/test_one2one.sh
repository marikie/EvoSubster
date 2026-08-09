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
printf '%s\n' "$@" > "$LASTDB_ARGS_LOG"
touch "${@: -2:1}.prj"
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

PATH="$stub_bin:$PATH" LASTAL_ARGS_LOG="$tmp_dir/lastal.args" \
    LASTDB_ARGS_LOG="$tmp_dir/lastdb.args" THREAD_NUM=4 \
    bash "$ROOT_DIR/src/align/one2one.sh" 20260809 "$reference" "$query" \
    "$tmp_dir/db" "$train" "$output" > "$tmp_dir/run.log" 2>&1

test -s "$output"
grep -Fxq -- '-i' "$tmp_dir/lastal.args"
grep -Fxq -- '1M' "$tmp_dir/lastal.args"
grep -Fxq -- '--split-b=8T' "$tmp_dir/lastal.args"
test "$(grep -Fxc -- '-P4' "$tmp_dir/lastal.args")" -eq 1
grep -Fxq -- '-uRY128' "$tmp_dir/lastdb.args"
test -d "$tmp_dir/db.seed-RY128"

rm -f "$output"
PATH="$stub_bin:$PATH" LASTAL_ARGS_LOG="$tmp_dir/lastal.args" \
    LASTDB_ARGS_LOG="$tmp_dir/lastdb.args" THREAD_NUM=2 \
    LASTDB_ALIGNMENT_SEED=RY4 LASTAL_SPLIT_MEMORY=32G \
    bash "$ROOT_DIR/src/align/one2one.sh" 20260809 "$reference" "$query" \
    "$tmp_dir/db" "$train" "$output" > "$tmp_dir/run-custom.log" 2>&1

grep -Fxq -- '--split-b=32G' "$tmp_dir/lastal.args"
grep -Fxq -- '-uRY4' "$tmp_dir/lastdb.args"
test -d "$tmp_dir/db.seed-RY4"

if PATH="$stub_bin:$PATH" LASTDB_ALIGNMENT_SEED='' \
    bash "$ROOT_DIR/src/align/one2one.sh" 20260809 "$reference" "$query" \
    "$tmp_dir/db" "$train" "$tmp_dir/empty-seed.maf" \
    > "$tmp_dir/run-empty-seed.log" 2>&1; then
    echo "one2one unexpectedly accepted an empty seed scheme" >&2
    exit 1
fi
grep -Fq 'LASTDB_ALIGNMENT_SEED must not be empty' "$tmp_dir/run-empty-seed.log"

cat > "$stub_bin/lastal" <<'EOF'
#!/bin/bash
echo 'last-split: skipping sequence chr1 (9999999999 bytes)' >&2
printf '##maf version=1\n\n'
EOF
chmod +x "$stub_bin/lastal"
rm -f "$output"
if PATH="$stub_bin:$PATH" THREAD_NUM=1 \
    bash "$ROOT_DIR/src/align/one2one.sh" 20260809 "$reference" "$query" \
    "$tmp_dir/db" "$train" "$output" > "$tmp_dir/run-skip.log" 2>&1; then
    echo "one2one unexpectedly accepted a skipped LAST sequence" >&2
    exit 1
fi
test ! -e "$output"
grep -Fq 'incomplete alignments are not accepted' "$tmp_dir/run-skip.log"

echo "ok: one2one uses a whole-genome seed, bounded query batches, and rejects skipped sequences"
