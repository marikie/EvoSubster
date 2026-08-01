#!/bin/bash

set -u

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/evosubster-wrapper-cache.XXXXXX")
trap 'rm -rf "$tmp_dir"' EXIT

fail=0
check() {
    local description="$1"
    shift
    if "$@"; then
        echo "ok: $description"
    else
        echo "FAIL: $description"
        fail=$((fail + 1))
    fi
}

# Build a minimal isolated checkout so absolute helper paths used by
# sbst_fromDwl.sh can be replaced with deterministic process-boundary stubs.
fixture_repo="$tmp_dir/repo"
stub_bin="$tmp_dir/bin"
mkdir -p "$fixture_repo/src/select" "$fixture_repo/src/lib" "$fixture_repo/config" "$stub_bin"
cp "$ROOT_DIR/src/sbst_fromDwl.sh" "$fixture_repo/src/sbst_fromDwl.sh"
cp "$ROOT_DIR/src/lib/train_cache.sh" "$fixture_repo/src/lib/train_cache.sh"
cp "$ROOT_DIR/config/dwl_config.yaml" "$fixture_repo/config/dwl_config.yaml"

cat > "$stub_bin/Rscript" <<'EOF'
#!/bin/bash
script="$1"
shift
printf '%s\n' "$script" "$@" > "$R_ARGS_LOG"
out_dir=""
while [ "$#" -gt 0 ]; do
    case "$1" in
        --out-dir)
            out_dir="$2"
            shift 2
            ;;
        *)
            shift
            ;;
    esac
done
mkdir -p "$out_dir"
printf 'out_acc\tin1_acc\tin2_acc\nGCA_000000011.1\tGCA_000000022.1\tGCA_000000033.1\n' \
    > "$out_dir/selected_trios.tsv"
EOF

cat > "$stub_bin/python" <<'EOF'
#!/bin/bash
exec python3 "$@"
EOF

cat > "$fixture_repo/src/dwl_organism.sh" <<'EOF'
#!/bin/bash
acc="$1"
shift
out_dir=""
while [ "$#" -gt 0 ]; do
    case "$1" in
        --out-dir)
            out_dir="$2"
            shift 2
            ;;
        *)
            shift
            ;;
    esac
done
if [[ "$acc" =~ \.[0-9]+$ ]]; then
    resolved_acc="$acc"
else
    resolved_acc="${acc}.9"
fi
org_name="Species_${resolved_acc//[._]/_}"
org_dir="$out_dir/$org_name"
mkdir -p "$org_dir"
fasta="$org_dir/${resolved_acc}_genomic.fna"
summary="$org_dir/summary.json"
taxonomy="$org_dir/taxonomy.json"
printf '>seq\nACGT\n' > "$fasta"
printf '{"reports":[{"accession":"%s"}]}\n' "$resolved_acc" > "$summary"
printf '{}\n' > "$taxonomy"
printf '%s|%s|%s|%s|%s|%s\n' \
    "$org_name" "$fasta" "$summary" "$org_name" "$org_name" "$taxonomy"
EOF

cat > "$fixture_repo/src/sbst.sh" <<'EOF'
#!/bin/bash
printf '%s\n' "$*" >> "$SBST_ARGS_LOG"
EOF

chmod +x \
    "$stub_bin/Rscript" \
    "$stub_bin/python" \
    "$fixture_repo/src/dwl_organism.sh" \
    "$fixture_repo/src/sbst.sh"

tree_file="$tmp_dir/tree.nwk"
printf '((A,B),C);\n' > "$tree_file"

run_tree_mode() {
    local output_dir="$1" custom_cache="${2:-}"
    local r_log="$3" sbst_log="$4"
    local args=(
        --tree "$tree_file"
        20260723
        --genome-dir "$tmp_dir/genomes"
        --out-dir "$output_dir"
        --idt-only
    )
    if [ -n "$custom_cache" ]; then
        args+=(--train-cache-dir "$custom_cache")
    fi

    PATH="$stub_bin:$PATH" \
    R_ARGS_LOG="$r_log" \
    SBST_ARGS_LOG="$sbst_log" \
    bash "$fixture_repo/src/sbst_fromDwl.sh" "${args[@]}" \
        > "$tmp_dir/wrapper.log" 2>&1
}

default_out="$tmp_dir/default-out"
default_r_log="$tmp_dir/default-r.log"
default_sbst_log="$tmp_dir/default-sbst.log"
run_tree_mode "$default_out" "" "$default_r_log" "$default_sbst_log"
default_exit=$?
default_cache="$default_out/trio_selection/train_cache"
if [ "$default_exit" -ne 0 ]; then
    sed 's/^/wrapper stderr: /' "$tmp_dir/wrapper.log"
fi

check "tree mode completes with the default cache path" test "$default_exit" -eq 0
check "default selector invocation relies on its documented cache default" \
    sh -c "! grep -q -- '--train-cache-dir' \"\$1\"" sh "$default_r_log"
check "default cache path is forwarded to the downstream pipeline" \
    grep -Fq -- "--train-cache-dir $default_cache" "$default_sbst_log"
check "selected accessions are forwarded explicitly to sbst.sh" \
    grep -Fq -- "--accession-ids GCA_000000011.1 GCA_000000022.1 GCA_000000033.1" \
    "$default_sbst_log"

custom_out="$tmp_dir/custom-out"
custom_cache="$tmp_dir/custom-cache"
custom_r_log="$tmp_dir/custom-r.log"
custom_sbst_log="$tmp_dir/custom-sbst.log"
run_tree_mode "$custom_out" "$custom_cache" "$custom_r_log" "$custom_sbst_log"
custom_exit=$?
if [ "$custom_exit" -ne 0 ]; then
    sed 's/^/wrapper stderr: /' "$tmp_dir/wrapper.log"
fi

check "tree mode completes with a custom cache path" test "$custom_exit" -eq 0
check "custom cache path is passed to the R selector" \
    sh -c 'grep -Fxq -- "$1" "$2"' sh "$custom_cache" "$custom_r_log"
check "the same custom cache path is passed to sbst.sh" \
    grep -Fq -- "--train-cache-dir $custom_cache" "$custom_sbst_log"

qc_file="$tmp_dir/external-qc.tsv"
printf 'accession\tqc_busco_mode\nGCA_000000011.1\tgenome\n' > "$qc_file"
quality_out="$tmp_dir/quality-out"
quality_r_log="$tmp_dir/quality-r.log"
quality_sbst_log="$tmp_dir/quality-sbst.log"
PATH="$stub_bin:$PATH" \
R_ARGS_LOG="$quality_r_log" \
SBST_ARGS_LOG="$quality_sbst_log" \
bash "$fixture_repo/src/sbst_fromDwl.sh" \
    --tree "$tree_file" 20260723 \
    --genome-dir "$tmp_dir/quality-genomes" \
    --out-dir "$quality_out" \
    --stage0-top-k 4 \
    --assembly-qc "$qc_file" \
    --idt-only \
    > "$tmp_dir/quality-wrapper.log" 2>&1
quality_exit=$?

check "tree wrapper accepts Stage 0 quality-ranking options" test "$quality_exit" -eq 0
check "tree wrapper forwards the Stage 0 shortlist size" \
    sh -c 'grep -Fxq -- "$1" "$2"' sh 4 "$quality_r_log"
check "tree wrapper forwards the external assembly-QC table" \
    sh -c 'grep -Fxq -- "$1" "$2"' sh "$qc_file" "$quality_r_log"

unversioned_log="$tmp_dir/unversioned-sbst.log"
PATH="$stub_bin:$PATH" \
SBST_ARGS_LOG="$unversioned_log" \
bash "$fixture_repo/src/sbst_fromDwl.sh" \
    20260723 GCA_000000011 GCA_000000022 GCA_000000033 \
    --genome-dir "$tmp_dir/unversioned-genomes" \
    --out-dir "$tmp_dir/unversioned-out" \
    --idt-only \
    > "$tmp_dir/unversioned-wrapper.log" 2>&1
unversioned_exit=$?

check "single-trio wrapper preserves support for unversioned requested accessions" \
    test "$unversioned_exit" -eq 0
check "single-trio wrapper forwards resolved versioned accessions" \
    grep -Fq -- "--accession-ids GCA_000000011.9 GCA_000000022.9 GCA_000000033.9" \
    "$unversioned_log"

if [ "$fail" -gt 0 ]; then
    echo
    echo "$fail test(s) FAILED"
    exit 1
fi

echo
echo "All sbst_fromDwl cache-forwarding tests passed."
