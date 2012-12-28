#!/bin/bash

# determine number of reads contained in fastq file (given as argument; possibly zipped)

# alternative is to work on FastQC data:
# awk  '/^Total Sequences/ {print $NF}' fastqc_data.txt
#
# http://en.wikipedia.org/wiki/FASTQ_format
#
# there are many ways to count the number of reads. in most cases you
# could also just divide the number of lines by four, however the spec
# does not require the sequence itself or the quality to be on one
# line zgrep -v '^@' might return hits on base quality lines as well
# and restricting it further is not possible since the fasta-id is
# optional. the following should never produce false positives as it
# also does some error checking.

test -z "$1" && exit 1
fastq="$1"
{ zcat $fastq 2>/dev/null || cat $fastq; } | \
  awk '{if ((NR-1)%4!=0) {next}; if (substr($0,0,1)!="@") {exit 1}; nreads+=1;} END {print nreads}'

