#!/bin/bash

# determine number of reads contained in fastq file (given as argument; possibly zipped)
#
# There are many ways to count the number of reads. in most cases you
# could also just divide the number of lines by four, however the spec
# does not require the sequence itself or the quality to be on just
# one line. furthermore zgrep -v '^@' might return hits on
# base-quality lines as well and restricting it further is not
# possible since the fasta-id is optional. the following also assumes
# that a read spans four lines, but it does at least perform some
# checks (and supports gzip)
#
# could check every fourth line in awk (if (NR%4 != 1) {next}) but
# that takes surprisingly longer than getting every fourth line with
# sed first (using non-gnu sed call for compatibility here. the
# following would work in gnu sed as well: -n '1~4p')

for fastq in "$@"; do
    echo -ne "$fastq\t"
    { zcat $fastq 2>/dev/null || cat $fastq; } | \
        #sed -n '1~4p' | \
        sed -n '1,${p;n;n;n;}' | \
        awk '{if (substr($0,0,1)!="@") {exit 1}; nreads+=1;} END {print nreads}'
done


