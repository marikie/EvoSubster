#!/usr/bin/env bash

set -u
set -o pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
SCRIPT="$ROOT_DIR/src/report/link_trio_dirs.sh"
TEST_ROOT=$(mktemp -d "${TMPDIR:-/tmp}/link-trio-dirs-test.XXXXXX")
trap 'rm -rf "$TEST_ROOT"' EXIT

failures=0

check() {
    local description="$1"
    shift
    if "$@"; then
        printf 'ok: %s\n' "$description"
    else
        printf 'FAIL: %s\n' "$description"
        failures=$((failures + 1))
    fi
}

make_trio() {
    local root="$1"
    local name="$2"
    mkdir -p "$root/$name/20260820"
    touch "$root/$name/20260820/marker"
}

run_linker() {
    bash "$SCRIPT" "$@"
}

lineage_a="$TEST_ROOT/actinopteri"
lineage_b="$TEST_ROOT/eupercaria"
output_dir="$TEST_ROOT/all"

make_trio "$lineage_a" "Fish1_Fish2_Fish3"
make_trio "$lineage_a" "Shared1_Shared2_Shared3"
mkdir -p "$lineage_a/pca" "$lineage_a/tmp" "$lineage_a/helper_A_B_C"
make_trio "$lineage_b" "Mammal1_Mammal2_Mammal3"
mkdir -p "$lineage_b/no_numeric_run_A_B_C"

first_log="$TEST_ROOT/first.log"
if run_linker --output-dir "$output_dir" "$lineage_a" "$lineage_b" >"$first_log" 2>&1; then
    first_status=0
else
    first_status=$?
fi
check "links multiple lineage roots successfully" test "$first_status" -eq 0
check "creates the first trio link" test -L "$output_dir/Fish1_Fish2_Fish3"
check "creates the second lineage trio link" test -L "$output_dir/Mammal1_Mammal2_Mammal3"
check "does not link pca or tmp directories" test ! -e "$output_dir/pca" -a ! -e "$output_dir/tmp"
check "does not link a named directory without a numeric run" test ! -e "$output_dir/helper_A_B_C" -a ! -e "$output_dir/no_numeric_run_A_B_C"
check "created link uses a relative target" test "$(readlink "$output_dir/Fish1_Fish2_Fish3")" = "../actinopteri/Fish1_Fish2_Fish3"
check "relative link resolves to the source directory" test "$(readlink -f "$output_dir/Fish1_Fish2_Fish3")" = "$(readlink -f "$lineage_a/Fish1_Fish2_Fish3")"
check "reports the number of created links" grep -q "Created links: 3" "$first_log"

second_log="$TEST_ROOT/second.log"
if run_linker --output-dir "$output_dir" "$lineage_a" "$lineage_b" >"$second_log" 2>&1; then
    second_status=0
else
    second_status=$?
fi
check "reuses correct existing links" test "$second_status" -eq 0
check "reports reused links" grep -q "Reused links: 3" "$second_log"

collision_a="$TEST_ROOT/collision_a"
collision_b="$TEST_ROOT/collision_b"
collision_output="$TEST_ROOT/collision_all"
make_trio "$collision_a" "Same1_Same2_Same3"
make_trio "$collision_b" "Same1_Same2_Same3"
collision_log="$TEST_ROOT/collision.log"
if run_linker --output-dir "$collision_output" "$collision_a" "$collision_b" >"$collision_log" 2>&1; then
    collision_status=0
else
    collision_status=$?
fi
check "fails on duplicate trio names" test "$collision_status" -ne 0
check "reports the duplicate trio" grep -q "collision" "$collision_log"
check "creates no links after a collision" test -z "$(find "$collision_output" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)"

wrong_output="$TEST_ROOT/wrong_link_all"
mkdir -p "$wrong_output"
ln -s "$lineage_b/Mammal1_Mammal2_Mammal3" "$wrong_output/Fish1_Fish2_Fish3"
wrong_log="$TEST_ROOT/wrong.log"
if run_linker --output-dir "$wrong_output" "$lineage_a" >"$wrong_log" 2>&1; then
    wrong_status=0
else
    wrong_status=$?
fi
check "fails on an existing link to the wrong target" test "$wrong_status" -ne 0

directory_output="$TEST_ROOT/directory_all"
mkdir -p "$directory_output/Fish1_Fish2_Fish3"
directory_log="$TEST_ROOT/directory.log"
if run_linker --output-dir "$directory_output" "$lineage_a" >"$directory_log" 2>&1; then
    directory_status=0
else
    directory_status=$?
fi
check "fails when the destination is a real directory" test "$directory_status" -ne 0

missing_log="$TEST_ROOT/missing.log"
if run_linker --output-dir "$TEST_ROOT/missing_output" "$TEST_ROOT/does-not-exist" >"$missing_log" 2>&1; then
    missing_status=0
else
    missing_status=$?
fi
check "fails for a missing input directory" test "$missing_status" -ne 0

if [ "$failures" -gt 0 ]; then
    exit 1
fi

printf 'All trio directory link tests passed.\n'
