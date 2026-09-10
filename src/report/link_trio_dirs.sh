#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'USAGE'
Usage:
  link_trio_dirs.sh --output-dir OUTPUT_DIR INPUT_DIR [INPUT_DIR...]

Collect eligible trio directories from each input lineage root and create
relative symlinks beneath OUTPUT_DIR.
USAGE
}

die() {
    printf 'Error: %s\n' "$*" >&2
    exit 1
}

has_numeric_run() {
    local trio_dir="$1"
    local child

    while IFS= read -r -d '' child; do
        if [[ "$(basename -- "$child")" =~ ^[0-9]+$ ]]; then
            return 0
        fi
    done < <(find "$trio_dir" -mindepth 1 -maxdepth 1 -type d -print0)

    return 1
}

is_trio_directory() {
    local candidate="$1"
    local name

    name=$(basename -- "$candidate")
    [[ "$name" =~ ^[^_]+_[^_]+_[^_]+$ ]] || return 1
    has_numeric_run "$candidate"
}

parse_args() {
    output_dir=''
    input_dirs=()

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --output-dir)
                [[ $# -ge 2 ]] || die "--output-dir requires a value"
                output_dir=$2
                shift 2
                ;;
            --output-dir=*)
                output_dir=${1#*=}
                shift
                ;;
            -h|--help)
                usage
                exit 0
                ;;
            --)
                shift
                while [[ $# -gt 0 ]]; do
                    input_dirs+=("$1")
                    shift
                done
                ;;
            -*)
                die "unknown option: $1"
                ;;
            *)
                input_dirs+=("$1")
                shift
                ;;
        esac
    done

    [[ -n "$output_dir" ]] || die "--output-dir is required"
    [[ ${#input_dirs[@]} -gt 0 ]] || die "at least one input directory is required"
}

parse_args "$@"

mkdir -p -- "$output_dir"
[[ -d "$output_dir" ]] || die "output directory is not a directory: $output_dir"
output_dir=$(cd -- "$output_dir" && pwd -P)

declare -A planned_sources=()
declare -a trio_names=()

for input_dir in "${input_dirs[@]}"; do
    [[ -d "$input_dir" ]] || die "input directory does not exist: $input_dir"
    input_dir=$(cd -- "$input_dir" && pwd -P)
    [[ "$input_dir" != "$output_dir" ]] || die "output directory cannot also be an input directory: $input_dir"

    while IFS= read -r -d '' candidate; do
        is_trio_directory "$candidate" || continue

        name=$(basename -- "$candidate")
        source_dir=$(cd -- "$candidate" && pwd -P)
        if [[ -v "planned_sources[$name]" ]]; then
            if [[ "${planned_sources[$name]}" != "$source_dir" ]]; then
                die "collision: trio '$name' exists in both '${planned_sources[$name]}' and '$source_dir'"
            fi
            continue
        fi

        planned_sources["$name"]=$source_dir
        trio_names+=("$name")
    done < <(find "$input_dir" -mindepth 1 -maxdepth 1 -type d -print0)
done

for name in "${trio_names[@]}"; do
    destination="$output_dir/$name"
    source_dir=${planned_sources[$name]}

    if [[ -L "$destination" ]]; then
        resolved_destination=$(readlink -f -- "$destination")
        [[ "$resolved_destination" == "$source_dir" ]] || die "existing link has the wrong target: $destination"
    elif [[ -e "$destination" ]]; then
        die "destination already exists and is not the expected symlink: $destination"
    fi
done

created_links=0
reused_links=0

for name in "${trio_names[@]}"; do
    destination="$output_dir/$name"
    source_dir=${planned_sources[$name]}

    if [[ -L "$destination" ]]; then
        reused_links=$((reused_links + 1))
        continue
    fi

    relative_target=$(realpath --relative-to="$output_dir" "$source_dir")
    ln -s -- "$relative_target" "$destination"
    created_links=$((created_links + 1))
done

printf 'Output directory: %s\n' "$output_dir"
printf 'Input directories: %d\n' "${#input_dirs[@]}"
printf 'Trio candidates: %d\n' "${#trio_names[@]}"
printf 'Created links: %d\n' "$created_links"
printf 'Reused links: %d\n' "$reused_links"
