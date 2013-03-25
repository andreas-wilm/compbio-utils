#!/bin/bash
awk '/^[^#]/ {s+=($3-$2)} END {print s}' $@
