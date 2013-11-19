#!/bin/bash

# give bed file as last arg or pipe in
awk '{if (/^[^#]/) {d[$1]+=($3-$2)}} END {for (c in d) {printf "%s\t%d\n", c, d[c]; s+=d[c]}; printf "TOTAL\t%d\n", s}' $@
