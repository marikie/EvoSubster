#!/bin/bash

set -u

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/evosubster-download.XXXXXX")
trap 'rm -rf "$tmp_dir"' EXIT

fixture_zip="$tmp_dir/ncbi_dataset.zip"
mkdir -p "$tmp_dir/archive/ncbi_dataset/data/GCA_000000001.1"
printf '>sequence\nACGT\n' \
    > "$tmp_dir/archive/ncbi_dataset/data/GCA_000000001.1/GCA_000000001.1_genomic.fna"
(
    cd "$tmp_dir/archive" || exit 1
    zip -qr "$fixture_zip" ncbi_dataset
) || exit 1

stub_bin="$tmp_dir/bin"
mkdir -p "$stub_bin"

cat > "$stub_bin/curl" <<'EOF'
#!/bin/bash
output=""
url=""
http1=0
while [ "$#" -gt 0 ]; do
    case "$1" in
        -o)
            output="$2"
            shift 2
            ;;
        --http1.1)
            http1=1
            shift
            ;;
        http*)
            url="$1"
            shift
            ;;
        *)
            shift
            ;;
    esac
done

if [ -n "$output" ]; then
    endpoint="genome"
elif [[ "$url" == *"/taxonomy/"* ]]; then
    endpoint="taxonomy"
else
    endpoint="summary"
fi

counter="$CURL_STATE_DIR/$endpoint"
attempt=1
if [ -f "$counter" ]; then
    attempt=$(( $(cat "$counter") + 1 ))
fi
printf '%s\n' "$attempt" > "$counter"
if [ "$http1" -eq 1 ]; then
    printf '%s\n' "$attempt" > "$CURL_STATE_DIR/$endpoint.http1"
fi
if [ "${CURL_ALWAYS_DNS:-0}" -eq 1 ] || [ "$attempt" -eq 1 ]; then
    exit "${CURL_FIRST_ERROR:-92}"
fi

if [ -n "$output" ]; then
    cp "$FIXTURE_ZIP" "$output"
elif [[ "$url" == *"/taxonomy/"* ]]; then
    printf '{"taxonomy":"test"}\n'
else
    printf '{"reports":[{"organism":{"organism_name":"Test species","tax_id":123}}]}\n'
fi
EOF

cat > "$stub_bin/sleep" <<'EOF'
#!/bin/bash
exit 0
EOF

cat > "$stub_bin/jq" <<'EOF'
#!/bin/bash
query="$*"
case "$query" in
    *infraspecific_names*) printf '\n' ;;
    *organism_name*) printf 'Test species\n' ;;
    *tax_id*) printf '123\n' ;;
    *) exit 1 ;;
esac
EOF

cat > "$stub_bin/yq" <<'EOF'
#!/bin/bash
case "$*" in
    *patterns.fasta*) printf '{org_id}*.fna\n' ;;
    *download.includes*) printf 'GENOME_FASTA\n' ;;
esac
EOF

chmod +x "$stub_bin/curl" "$stub_bin/jq" "$stub_bin/sleep" "$stub_bin/yq"

work_dir="$tmp_dir/work"
curl_state_dir="$tmp_dir/curl-state"
mkdir -p "$work_dir" "$curl_state_dir"
result=$(
    cd "$work_dir" || exit 1
    PATH="$stub_bin:$PATH" FIXTURE_ZIP="$fixture_zip" CURL_STATE_DIR="$curl_state_dir" \
        bash "$ROOT_DIR/src/dwl_organism.sh" GCA_000000001.1 --out-dir ./genomes
) || exit 1

IFS='|' read -r dir_name fasta_path summary_path raw_name full_name taxonomy_path <<< "$result"

test "$dir_name" = "Test_species"
test "$fasta_path" = "./genomes/Test_species/GCA_000000001.1_genomic.fna"
test -s "$work_dir/$fasta_path"
test -s "$work_dir/$summary_path"
test -s "$work_dir/$taxonomy_path"
test "$raw_name" = "Test species"
test "$full_name" = "Test_species"
test "$(cat "$curl_state_dir/summary")" -eq 2
test "$(cat "$curl_state_dir/genome")" -eq 2
test "$(cat "$curl_state_dir/taxonomy")" -eq 2
test "$(cat "$curl_state_dir/summary.http1")" -eq 2
test "$(cat "$curl_state_dir/genome.http1")" -eq 2
test "$(cat "$curl_state_dir/taxonomy.http1")" -eq 2

failure_state_dir="$tmp_dir/failure-state"
mkdir -p "$failure_state_dir"
if (
    cd "$work_dir" || exit 1
    PATH="$stub_bin:$PATH" FIXTURE_ZIP="$fixture_zip" CURL_STATE_DIR="$failure_state_dir" \
        CURL_ALWAYS_DNS=1 bash "$ROOT_DIR/src/dwl_organism.sh" GCF_000000002.1 \
        --out-dir ./genomes --no-genome --no-taxonomy \
        > "$tmp_dir/failure.stdout" 2> "$tmp_dir/failure.stderr"
); then
    echo "FAIL: sustained DNS failure unexpectedly succeeded" >&2
    exit 1
fi
test "$(cat "$failure_state_dir/summary")" -eq 4
grep -Fq "Error: Failed to fetch genome summary for GCF_000000002.1" \
    "$tmp_dir/failure.stderr"

echo "ok: fresh download retries HTTP/2 failures over HTTP/1.1 and preserves relative output paths"
echo "ok: sustained DNS failure stops after four attempts"