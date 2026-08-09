#!/bin/bash

set -eu

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/evosubster-atomic.XXXXXX")
trap 'rm -rf "$tmp_dir"' EXIT
stub_bin="$tmp_dir/bin"
mkdir -p "$stub_bin"

cat > "$stub_bin/maf-sort" <<'EOF'
#!/bin/bash
cat "$1"
EOF
cat > "$stub_bin/maf-join" <<'EOF'
#!/bin/bash
printf 'partial\n'
[ "${FAIL_MAF_JOIN:-0}" -eq 0 ]
EOF
chmod +x "$stub_bin/maf-sort" "$stub_bin/maf-join"

printf 'maf-a\n' > "$tmp_dir/a.maf"
printf 'maf-b\n' > "$tmp_dir/b.maf"
if PATH="$stub_bin:$PATH" FAIL_MAF_JOIN=1 \
    bash "$ROOT_DIR/src/align/mafjoin.sh" \
    "$tmp_dir/a.maf" "$tmp_dir/b.maf" "$tmp_dir/joined.maf" >/dev/null 2>&1; then
    echo "FAIL: interrupted maf-join unexpectedly succeeded" >&2
    exit 1
fi
test ! -e "$tmp_dir/joined.maf"

PATH="$stub_bin:$PATH" bash "$ROOT_DIR/src/align/mafjoin.sh" \
    "$tmp_dir/a.maf" "$tmp_dir/b.maf" "$tmp_dir/joined.maf" >/dev/null 2>&1
test -s "$tmp_dir/joined.maf"
test -e "$tmp_dir/joined.maf.complete"

cat > "$stub_bin/python3" <<'EOF'
#!/bin/bash
while [ "$#" -gt 0 ]; do
    case "$1" in
        -o2) out2="$2"; shift 2 ;;
        -o3) out3="$2"; shift 2 ;;
        *) shift ;;
    esac
done
printf 'partial\n' > "$out2"
printf 'partial\n' > "$out3"
[ "${FAIL_TSV:-0}" -eq 0 ]
EOF
chmod +x "$stub_bin/python3"

if PATH="$stub_bin:$PATH" FAIL_TSV=1 \
    bash "$ROOT_DIR/src/generate_tsv_files.sh" \
    "$tmp_dir/joined.maf" "$tmp_dir/out2.tsv" "$tmp_dir/out3.tsv" "$tmp_dir" \
    "$tmp_dir/out2-dinuc.tsv" "$tmp_dir/out3-dinuc.tsv" >/dev/null 2>&1; then
    echo "FAIL: interrupted TSV generation unexpectedly succeeded" >&2
    exit 1
fi
test ! -e "$tmp_dir/out2.tsv"
test ! -e "$tmp_dir/out3.tsv"

PATH="$stub_bin:$PATH" bash "$ROOT_DIR/src/generate_tsv_files.sh" \
    "$tmp_dir/joined.maf" "$tmp_dir/out2.tsv" "$tmp_dir/out3.tsv" "$tmp_dir" \
    "$tmp_dir/out2-dinuc.tsv" "$tmp_dir/out3-dinuc.tsv" >/dev/null 2>&1
test -s "$tmp_dir/out2.tsv"
test -s "$tmp_dir/out3.tsv"
test -e "$tmp_dir/out2.tsv.complete"

echo "ok: interrupted MAF joins and TSV generation leave no reusable partial output"
