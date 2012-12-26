#!/bin/bash
#
# check wether quality encoding in fastq file (arg) is illumina or
# sanger. extracts quality encoding from pre-run fastqc.
fastq="$1"
test -n "$fastq" || exit 1

# infer fastqc data filename: supported fastq extensions are fastq and
# txt with or without gzip
fastqc_data=$(ls $(echo ${fastq%.gz} | \
    sed -e 's,\(.fastq\|.txt\),,')_fastqc/fastqc_data.txt) || exit 1
awk '/^Encoding/ {
  if ($NF>=1.8 || /Sanger/) {print "sanger"} else {print "illumina"}
}' $fastqc_data 
