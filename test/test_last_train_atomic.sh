#!/bin/bash

set -eu

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/evosubster-last-train-atomic.XXXXXX")
trap 'rm -rf "$tmp_dir"' EXIT

stub_bin="$tmp_dir/bin"
mkdir -p "$stub_bin"

reference="$tmp_dir/reference.fna"
query="$tmp_dir/query.fna"
printf '>reference\nACGT\n' > "$reference"
printf '>query\nACGT\n' > "$query"

cat > "$stub_bin/lastdb" <<'EOF'
#!/bin/bash
touch partial.prj
exit 1
EOF

cat > "$stub_bin/last-train" <<'EOF'
#!/bin/bash
printf '# partial parameters\n'
exit 1
EOF

chmod +x "$stub_bin/lastdb" "$stub_bin/last-train"

failed_db="$tmp_dir/failed-db"
if PATH="$stub_bin:$PATH" THREAD_NUM=1 \
    bash "$ROOT_DIR/src/align/last_train.sh" 20260809 "$reference" "$query" Ref Qry "$failed_db" \
    > "$tmp_dir/failed-db.log" 2>&1; then
    echo "last_train unexpectedly accepted a failed lastdb build" >&2
    exit 1
fi
test ! -e "$failed_db/Refdb_20260809"
test -z "$(find "$failed_db" -maxdepth 1 -name 'Refdb_20260809.tmp.*' -print -quit)"

cat > "$stub_bin/lastdb" <<'EOF'
#!/bin/bash
printf 'complete\n' > "${@: -2:1}.prj"
EOF
chmod +x "$stub_bin/lastdb"

failed_train="$tmp_dir/failed-train"
if PATH="$stub_bin:$PATH" THREAD_NUM=1 \
    bash "$ROOT_DIR/src/align/last_train.sh" 20260809 "$reference" "$query" Ref Qry "$failed_train" \
    > "$tmp_dir/failed-train.log" 2>&1; then
    echo "last_train unexpectedly accepted a failed parameter fit" >&2
    exit 1
fi
test -s "$failed_train/Refdb_20260809/Refdb_20260809.prj"
test ! -e "$failed_train/Ref2Qry_20260809.train"
test -z "$(find "$failed_train" -maxdepth 1 -name 'Ref2Qry_20260809.train.tmp.*' -print -quit)"

cat > "$stub_bin/last-train" <<'EOF'
#!/bin/bash
printf '# complete parameters\n'
EOF
chmod +x "$stub_bin/last-train"

success_dir="$tmp_dir/success"
PATH="$stub_bin:$PATH" THREAD_NUM=1 \
    bash "$ROOT_DIR/src/align/last_train.sh" 20260809 "$reference" "$query" Ref Qry "$success_dir" \
    > "$tmp_dir/success.log" 2>&1

test -s "$success_dir/Refdb_20260809/Refdb_20260809.prj"
test -s "$success_dir/Ref2Qry_20260809.train"

echo "ok: last_train publishes databases and parameter files only after successful completion"
