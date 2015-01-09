#!/bin/bash

bam="$1"
test -n "$bam" || exit 1
samtools view -H $bam | \
	awk '/^@SQ/ {
if (substr($2, 1, 3)!="SN:" || substr($3, 0, 3)!="LN:") {exit 1;}
sq=substr($2, 4, length($2)-3); ln=substr($3, 4, length($3)-3); printf "%s\t%d\t%d\n", sq, 0, ln}' 
