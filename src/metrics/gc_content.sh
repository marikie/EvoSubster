#!/bin/bash

numer=$(cat "$@" | grep -v '>' | tr -cd CGcg | wc -c)
denom=$(cat "$@" | grep -v '>' | tr -cd ACGTacgt | wc -c)

awk 'BEGIN {printf "Total GC content: %.2f%%\n", '$numer' / '$denom' * 100}'