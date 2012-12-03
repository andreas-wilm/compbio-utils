#!/bin/bash

# extract barcode from fastq file (gzip supported).  default is to check
# only the first line, i.e.  only one barcode will be extracted.  an
# optional 2nd argument defining the max number of reads to check, can be
# given.
#
# assuming each read spans four lines

fq="$1"
test -n "$fq" || exit 1
test -s "$fq" || exit 1

# optional 2nd arg: max number of lines to check
ml=1
test -s "$ml" || ml=$2

# use cat or zcat?
cat=cat
file "$fq" | grep -q gzip && cat=zcat

# element after last hash or colon is barcode, possibly with mate pair. 
# should work with all fastq files. note, using sort -u instead of awk
# array for sake of readability.
$cat "$fq" | \
    awk -v ml=$ml '{if (NR%4==1) {sub(/.*[#:]/, ""); sub(/\/[12]$/, ""); print; s+=1; if (s>ml) {exit}}}' | sort -u

