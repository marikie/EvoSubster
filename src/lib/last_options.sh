#!/bin/bash

LAST_OPTIONS_HELPER="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/last_options.py"

last_options_parse_cli() {
    LAST_REMAINING_ARGS=()
    local name value json
    json="${EVOSUBSTER_LAST_OPTIONS:-}"
    [ -n "$json" ] || json='{}'
    local -A seen=()
    while [ "$#" -gt 0 ]; do
        case "$1" in
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
                json=$(python3 -c 'import json,sys; d=json.loads(sys.argv[1]); d[sys.argv[2]]=sys.argv[3]; print(json.dumps(d))' "$json" "$name" "$value") || return 1
                ;;
            --)
                LAST_REMAINING_ARGS+=("$@")
                break
                ;;
            *) LAST_REMAINING_ARGS+=("$1"); shift ;;
        esac
    done
    export EVOSUBSTER_LAST_OPTIONS="$json"
    # Validate before downloading genomes or starting tree selection.
    python3 "$LAST_OPTIONS_HELPER" validate >/dev/null
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
