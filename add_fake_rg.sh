#!/bin/bash

# Add minimal fake default read group to a BAM file. GATK requires
# read groups (you wonder why there's no option to assume a default
# one, but never mind).

bam="$1";
test -n "$bam" || exit 1
test -s "$bam" || exit 1

samplename="$2"
test -n "$samplename" || samplename="sample-name-dummy"

# lazy implementation: parse twice. First head and add RG declaration.
# Then assign RG to reads.
samtools view -H $bam || exit 1
echo -e "@RG\tID:1\tPL:illumina\tPU:platform-unit-dummy\tLB:library-dummy\tSM:${samplename}\tCN:sequencing-center-dummy";
samtools view $bam | \
     awk 'BEGIN {FS="\t"; OFS="\t"} {printf "%s\tRG:Z:1\n", $0}';
