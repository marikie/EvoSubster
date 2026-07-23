#!/bin/bash

extract_ncbi_accession_from_path() {
    local base
    base=$(basename "$1")
    if [[ "$base" =~ ^((GCA|GCF)_[0-9]+\.[0-9]+)(_|\.|$) ]]; then
        printf '%s\n' "${BASH_REMATCH[1]}"
    fi
}

is_valid_train_file() {
    local path="$1"
    [ -f "$path" ] && [ -s "$path" ] || return 1

    # last-train prints several intermediate identity/score blocks. Only the
    # final block contains the #last parameters and uncommented 4x4 matrix that
    # lastal -p needs. Reset on every identity line so a truncated intermediate
    # block can never be mistaken for a complete parameter file.
    awk '
        function is_identity_number(value) {
            return value ~ /^[0-9]+([.][0-9]+)?$/ &&
                   (value + 0) >= 0 && (value + 0) <= 100
        }
        function is_score_number(value) {
            return value ~ /^-?[0-9]+([.][0-9]+)?$/
        }
        /^# substitution percent identity:[[:space:]]*/ {
            value = $0
            sub(/^# substitution percent identity:[[:space:]]*/, "", value)
            valid_identity = is_identity_number(value)
            seen_identity = 1
            have_t = have_a = have_A = have_b = have_B = have_S = 0
            in_matrix = 0
            delete matrix_rows
            next
        }
        !seen_identity { next }
        /^#last -t/ { have_t = 1; next }
        /^#last -a/ { have_a = 1; next }
        /^#last -A/ { have_A = 1; next }
        /^#last -b/ { have_b = 1; next }
        /^#last -B/ { have_B = 1; next }
        /^#last -S/ { have_S = 1; next }
        /^# score matrix \(query letters = columns, reference letters = rows\):$/ {
            in_matrix = 1
            next
        }
        in_matrix && $1 ~ /^[ACGT]$/ && NF == 5 {
            if (is_score_number($2) && is_score_number($3) &&
                is_score_number($4) && is_score_number($5)) {
                matrix_rows[$1] = 1
            }
        }
        END {
            complete = valid_identity && have_t && have_a && have_A &&
                       have_b && have_B && have_S &&
                       ("A" in matrix_rows) && ("C" in matrix_rows) &&
                       ("G" in matrix_rows) && ("T" in matrix_rows)
            exit !complete
        }
    ' "$path"
}

copy_train_file_atomically() {
    local source="$1" dest="$2"
    local tmp="${dest}.cache-copy.$$"

    rm -f "$tmp"
    if ! cp "$source" "$tmp"; then
        rm -f "$tmp"
        return 1
    fi
    if ! is_valid_train_file "$tmp"; then
        echo "Warning: copied train cache failed validation: $source" >&2
        rm -f "$tmp"
        return 1
    fi
    if ! mv -f "$tmp" "$dest"; then
        rm -f "$tmp"
        return 1
    fi
}

reuse_cached_train() {
    local cache_dir="$1" date="$2" ref_acc="$3" query_acc="$4" dest="$5"

    # Materialize valid symlinks from older runs so later cache cleanup cannot
    # break an otherwise complete pipeline output.
    if [ -L "$dest" ] && is_valid_train_file "$dest"; then
        if copy_train_file_atomically "$dest" "$dest"; then
            echo "Materialized legacy train-cache symlink: $dest"
            return 0
        fi
        echo "Error: failed to materialize legacy train-cache symlink: $dest" >&2
        return 1
    fi

    # last_train.sh skips every existing path, including partial output. Remove
    # invalid regular files and dangling symlinks before cache lookup so a fresh
    # training run can recover safely.
    if [ -e "$dest" ] || [ -L "$dest" ]; then
        if is_valid_train_file "$dest"; then
            return 0
        fi
        echo "Removing invalid train destination before retry: $dest" >&2
        if ! rm -f "$dest"; then
            echo "Error: failed to remove invalid train destination: $dest" >&2
            return 1
        fi
    fi

    if [ -z "$cache_dir" ] || [ -z "$ref_acc" ] || [ -z "$query_acc" ]; then
        return 0
    fi

    if [[ "$query_acc" < "$ref_acc" ]]; then
        echo "Train-cache skip: $ref_acc/$query_acc direction does not match the selector cache convention"
        return 0
    fi

    local cached="${cache_dir}/${ref_acc}2${query_acc}_${date}.train"
    if ! is_valid_train_file "$cached"; then
        if [ -e "$cached" ] || [ -L "$cached" ]; then
            echo "Train-cache miss: invalid or incomplete file at $cached" >&2
        else
            echo "Train-cache miss: no cached file at $cached"
        fi
        return 0
    fi

    if copy_train_file_atomically "$cached" "$dest"; then
        echo "Reused cached train file: $cached -> $dest"
        return 0
    fi

    # A failed copy must not leave a partial destination that would suppress
    # last-train. Falling back is safe because copy_train_file_atomically writes
    # to a temporary path first.
    echo "Warning: failed to copy cached train file; running last-train instead: $cached" >&2
    return 0
}
