#!/bin/bash

set -u

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TRAIN_CACHE_LIB="$ROOT_DIR/src/lib/train_cache.sh"

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

if [ ! -f "$TRAIN_CACHE_LIB" ]; then
    echo "FAIL: train-cache helper exists at $TRAIN_CACHE_LIB"
    exit 1
fi
# shellcheck source=../src/lib/train_cache.sh
source "$TRAIN_CACHE_LIB"

tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/evosubster-train-cache.XXXXXX")
trap 'rm -rf "$tmp_dir"' EXIT

write_valid_train() {
    local path="$1" identity="${2:-95.5}"
    printf '# substitution percent identity: %s\nscore matrix\n' "$identity" > "$path"
}

# Strict accession extraction must never turn an arbitrary FASTA name into a cache key.
check "extracts a complete NCBI accession from a downloader-style filename" \
    test "$(extract_ncbi_accession_from_path "/genomes/GCA_012345678.9_genomic.fna")" = \
    "GCA_012345678.9"
check "rejects a custom FASTA filename as an accession source" \
    test -z "$(extract_ncbi_accession_from_path "/genomes/custom_species.fna")"
check "rejects an accession without a version" \
    test -z "$(extract_ncbi_accession_from_path "/genomes/GCF_012345678_genomic.fna")"

cache_dir="$tmp_dir/cache"
dest_dir="$tmp_dir/dest"
mkdir -p "$cache_dir" "$dest_dir"

# Forward-orientation reuse copies validated content and survives cache cleanup.
forward_cache="$cache_dir/GCA_000000001.12GCA_000000002.1_20260723.train"
forward_dest="$dest_dir/forward.train"
write_valid_train "$forward_cache"
reuse_cached_train "$cache_dir" "20260723" \
    "GCA_000000001.1" "GCA_000000002.1" "$forward_dest"
check "reuses a forward-orientation cache entry" test -f "$forward_dest"
check "materializes a regular file instead of a symlink" test ! -L "$forward_dest"
rm -f "$forward_cache"
check "reused output remains valid after cache cleanup" is_valid_train_file "$forward_dest"

# Reverse orientation must not use parameters trained in the opposite direction.
reverse_cache="$cache_dir/GCA_000000001.12GCA_000000002.1_20260723.train"
reverse_dest="$dest_dir/reverse.train"
write_valid_train "$reverse_cache"
reuse_cached_train "$cache_dir" "20260723" \
    "GCA_000000002.1" "GCA_000000001.1" "$reverse_dest"
check "skips a reverse-orientation cache entry" test ! -e "$reverse_dest"

# A dangling destination from an older symlink-based run must be recoverable.
dangling_dest="$dest_dir/dangling.train"
ln -s "$tmp_dir/missing-cache.train" "$dangling_dest"
reuse_cached_train "$cache_dir" "20260723" \
    "GCA_000000001.1" "GCA_000000002.1" "$dangling_dest"
check "replaces a dangling destination link" test -f "$dangling_dest"
check "recovered destination is no longer a symlink" test ! -L "$dangling_dest"

# A valid symlink created by the previous implementation is converted in place
# so upgrading the code also makes existing outputs independent of the cache.
legacy_cache="$cache_dir/GCA_000000005.12GCA_000000006.1_20260723.train"
legacy_dest="$dest_dir/legacy.train"
write_valid_train "$legacy_cache"
ln -s "$legacy_cache" "$legacy_dest"
reuse_cached_train "$cache_dir" "20260723" \
    "GCA_000000005.1" "GCA_000000006.1" "$legacy_dest"
rm -f "$legacy_cache"
check "materializes a valid legacy cache symlink" test ! -L "$legacy_dest"
check "materialized legacy output survives cache cleanup" is_valid_train_file "$legacy_dest"

# Partial cache files must not be reused or leave a destination that makes last_train skip.
invalid_cache="$cache_dir/GCA_000000003.12GCA_000000004.1_20260723.train"
invalid_dest="$dest_dir/invalid.train"
printf 'partial output\n' > "$invalid_cache"
printf 'partial destination\n' > "$invalid_dest"
reuse_cached_train "$cache_dir" "20260723" \
    "GCA_000000003.1" "GCA_000000004.1" "$invalid_dest"
check "removes an invalid destination before fresh training" test ! -e "$invalid_dest"

# sbst.sh must use accessions supplied by its downloader caller even when FASTA
# basenames are custom, and it must attempt cache reuse for all three pairs.
stub_bin="$tmp_dir/bin"
stub_align="$tmp_dir/align"
stub_metrics="$tmp_dir/metrics"
pipeline_out="$tmp_dir/pipeline-out"
pipeline_cache="$tmp_dir/pipeline-cache"
pipeline_log="$tmp_dir/last-train.log"
mkdir -p "$stub_bin" "$stub_align" "$stub_metrics" "$pipeline_out" "$pipeline_cache"

cat > "$stub_bin/lastal" <<'EOF'
#!/bin/bash
echo "lastal test stub"
EOF
cat > "$stub_bin/readlink" <<'EOF'
#!/bin/bash
echo "readlink -f is unavailable in this compatibility test" >&2
exit 64
EOF
cat > "$stub_metrics/gc_content.sh" <<'EOF'
#!/bin/bash
echo "GC 0.5"
EOF
cat > "$stub_align/last_train.sh" <<'EOF'
#!/bin/bash
dest="$6/${4}2${5}_${1}.train"
if [ -s "$dest" ] && grep -Fq '# substitution percent identity:' "$dest"; then
    echo "reused" >> "$TRAIN_TEST_LOG"
else
    echo "fresh" >> "$TRAIN_TEST_LOG"
    printf '# substitution percent identity: 90\n' > "$dest"
fi
EOF
chmod +x \
    "$stub_bin/lastal" \
    "$stub_bin/readlink" \
    "$stub_metrics/gc_content.sh" \
    "$stub_align/last_train.sh"

mkdir -p "$tmp_dir/Genus_alpha" "$tmp_dir/Genus_beta" "$tmp_dir/Genus_gamma"
fasta1="$tmp_dir/Genus_alpha/custom-outgroup.fa"
fasta2="$tmp_dir/Genus_beta/custom-ingroup-a.fa"
fasta3="$tmp_dir/Genus_gamma/custom-ingroup-b.fa"
printf '>a\nACGT\n' > "$fasta1"
printf '>b\nACGT\n' > "$fasta2"
printf '>c\nACGT\n' > "$fasta3"

acc1="GCA_000000011.1"
acc2="GCA_000000022.1"
acc3="GCA_000000033.1"
write_valid_train "$pipeline_cache/${acc1}2${acc2}_20260723.train" 91
write_valid_train "$pipeline_cache/${acc1}2${acc3}_20260723.train" 92
write_valid_train "$pipeline_cache/${acc2}2${acc3}_20260723.train" 93

PATH="$stub_bin:$PATH" \
ALIGN_DIR_OVERRIDE="$stub_align" \
METRICS_DIR_OVERRIDE="$stub_metrics" \
TRAIN_TEST_LOG="$pipeline_log" \
bash "$ROOT_DIR/src/sbst.sh" \
    --train-cache-dir "$pipeline_cache" \
    --accession-ids "$acc1" "$acc2" "$acc3" \
    --out-dir "$pipeline_out" \
    --idt-only \
    20260723 "$fasta1" "$fasta2" "$fasta3" NO_GFF_FILE \
    > "$tmp_dir/sbst.log" 2>&1
pipeline_status=$?

check "sbst.sh accepts explicit accessions for custom FASTA names" \
    test "$pipeline_status" -eq 0
reuse_count=$(grep -c '^reused$' "$pipeline_log" 2>/dev/null || true)
fresh_count=$(grep -c '^fresh$' "$pipeline_log" 2>/dev/null || true)
reuse_count=${reuse_count:-0}
fresh_count=${fresh_count:-0}
check "sbst.sh reuses cached parameters for all three pairwise trainings" \
    test "$reuse_count" -eq 3
check "cache reuse does not require readlink -f" \
    sh -c "! grep -q 'readlink -f is unavailable' \"\$1\"" sh "$tmp_dir/sbst.log"
check "explicit-accession cache hits avoid every fresh last-train run" \
    test "$pipeline_status" -eq 0 -a "$fresh_count" -eq 0

if [ "$fail" -gt 0 ]; then
    echo
    echo "$fail test(s) FAILED"
    exit 1
fi

echo
echo "All train-cache tests passed."
