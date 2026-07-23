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
    [ -f "$path" ] &&
        [ -s "$path" ] &&
        grep -Fq '# substitution percent identity:' "$path"
}

copy_train_file_atomically() {
    local source="$1" dest="$2"
    local tmp="${dest}.cache-copy.$$"

    rm -f "$tmp"
    if ! cp "$source" "$tmp"; then
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
