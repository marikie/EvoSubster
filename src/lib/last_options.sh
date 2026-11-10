#!/bin/bash

LAST_OPTIONS_HELPER="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/last_options.py"

last_options_parse_cli() {
    LAST_REMAINING_ARGS=()
    LAST_EXPLICIT_NAMES=()
    local name value json
    local lastal_threads="${LASTAL_THREADS-}" threads_set="${LASTAL_THREADS+x}"
    local alignment_mode="${LASTAL_PIPELINE_MODE-streaming}"
    local reuse_from="${EVOSUBSTER_REUSE_FROM-}"
    local reuse_last_version="${EVOSUBSTER_REUSE_LAST_VERSION-same}"
    local reuse_plan="${EVOSUBSTER_REUSE_PLAN-0}"
    json="${EVOSUBSTER_LAST_OPTIONS:-}"
    [ -n "$json" ] || json='{}'
    local -A seen=()
    while [ "$#" -gt 0 ]; do
        case "$1" in
            --lastal-threads|--alignment-mode|--reuse-from|--reuse-last-version|\
            --lastal-threads=*|--alignment-mode=*|--reuse-from=*|--reuse-last-version=*)
                name="${1%%=*}"
                if [ -n "${seen[$name]:-}" ]; then
                    echo "Error: $name may only be specified once." >&2
                    return 1
                fi
                seen[$name]=1
                if [[ "$1" == *=* ]]; then
                    value="${1#*=}"
                    shift
                else
                    if [ "$#" -lt 2 ] || [[ "$2" == --* ]]; then
                        echo "Error: $name requires a value." >&2
                        return 1
                    fi
                    value="$2"
                    shift 2
                fi
                if [ -z "$value" ]; then
                    echo "Error: $name requires a non-empty value." >&2
                    return 1
                fi
                case "$name" in
                    --lastal-threads) lastal_threads="$value"; threads_set=x ;;
                    --alignment-mode) alignment_mode="$value" ;;
                    --reuse-from) reuse_from="$value" ;;
                    --reuse-last-version) reuse_last_version="$value" ;;
                esac
                ;;
            --reuse-plan)
                if [ -n "${seen[--reuse-plan]:-}" ]; then
                    echo "Error: --reuse-plan may only be specified once." >&2
                    return 1
                fi
                seen[--reuse-plan]=1
                reuse_plan=1
                shift
                ;;
            --reuse-plan=*)
                echo "Error: --reuse-plan does not take a value." >&2
                return 1
                ;;
            --lastdb-options|--lastal-options|--last-train-options|--last-split-options|\
            --lastdb-options=*|--lastal-options=*|--last-train-options=*|--last-split-options=*)
                name="${1%%=*}"
                if [ -n "${seen[$name]:-}" ]; then
                    echo "Error: $name may only be specified once." >&2
                    return 1
                fi
                seen[$name]=1
                if [[ "$1" == *=* ]]; then
                    value="${1#*=}"
                    shift
                else
                    if [ "$#" -lt 2 ]; then
                        echo "Error: $name requires an option string." >&2
                        return 1
                    fi
                    value="$2"
                    shift 2
                fi
                name="${name#--}"
                name="${name%-options}"
                LAST_EXPLICIT_NAMES+=("$name")
                json=$(python3 -c 'import json,sys; d=json.loads(sys.argv[1]); d[sys.argv[2]]=sys.argv[3]; print(json.dumps(d))' "$json" "$name" "$value") || return 1
                ;;
            --)
                LAST_REMAINING_ARGS+=("$@")
                break
                ;;
            *) LAST_REMAINING_ARGS+=("$1"); shift ;;
        esac
    done
    # Validate public runtime controls even when LAST option validation is deferred
    # until the entrypoint has resolved --thread. Do not create source directories.
    if [ -n "$threads_set" ] && ! [[ "$lastal_threads" =~ ^[1-9][0-9]*$ ]]; then
        echo "Error: --lastal-threads (LASTAL_THREADS) must be a positive integer." >&2
        return 1
    fi
    case "$alignment_mode" in
        streaming|staged) ;;
        *) echo "Error: --alignment-mode (LASTAL_PIPELINE_MODE) must be streaming or staged." >&2; return 1 ;;
    esac
    case "$reuse_last_version" in
        same|allow-existing-pairs) ;;
        *) echo "Error: --reuse-last-version (EVOSUBSTER_REUSE_LAST_VERSION) must be same or allow-existing-pairs." >&2; return 1 ;;
    esac
    case "$reuse_plan" in
        0|1) ;;
        *) echo "Error: --reuse-plan (EVOSUBSTER_REUSE_PLAN) must be 0 or 1." >&2; return 1 ;;
    esac
    if [ -n "$reuse_from" ]; then
        if [ ! -d "$reuse_from" ]; then
            echo "Error: --reuse-from (EVOSUBSTER_REUSE_FROM) must name an existing directory: $reuse_from" >&2
            return 1
        fi
        reuse_from=$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$reuse_from") || return 1
    fi
    [ -z "$threads_set" ] || export LASTAL_THREADS="$lastal_threads"
    export LASTAL_PIPELINE_MODE="$alignment_mode"
    export EVOSUBSTER_REUSE_FROM="$reuse_from"
    export EVOSUBSTER_REUSE_LAST_VERSION="$reuse_last_version"
    export EVOSUBSTER_REUSE_PLAN="$reuse_plan"
    export EVOSUBSTER_LAST_OPTIONS="$json"
    # Validate before downloading genomes or starting tree selection.
    if [ "${LAST_OPTIONS_DEFER_VALIDATE:-0}" != 1 ]; then
        python3 "$LAST_OPTIONS_HELPER" validate >/dev/null
    fi
}

last_options_load() {
    python3 "$LAST_OPTIONS_HELPER" validate >/dev/null || return 1
    mapfile -d '' -t LASTDB_TRAIN_ARGS < <(python3 "$LAST_OPTIONS_HELPER" argv lastdb_train)
    mapfile -d '' -t LASTDB_ALIGN_ARGS < <(python3 "$LAST_OPTIONS_HELPER" argv lastdb_align)
    mapfile -d '' -t LAST_TRAIN_ARGS < <(python3 "$LAST_OPTIONS_HELPER" argv last_train)
    mapfile -d '' -t LASTAL_ARGS < <(python3 "$LAST_OPTIONS_HELPER" argv lastal)
    mapfile -d '' -t LAST_SPLIT_ARGS < <(python3 "$LAST_OPTIONS_HELPER" argv last_split)
}

last_options_guard() {
    python3 "$LAST_OPTIONS_HELPER" guard "$1" "$2"
}
