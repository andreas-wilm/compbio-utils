#!/bin/bash
# create a "genome file" for genomeCoverageBed (why does that thing not support bed?)

bam=$1
test -z "$bam" && exit 1
test -e "$bam" || exit 1
samtools view -H $bam | grep '^@SQ' | awk '{printf "%s\t%s\n", substr($2, 4), substr($3, 4)}'

