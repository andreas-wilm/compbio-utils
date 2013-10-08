#!/bin/bash

# Author: Andreas Wilm <andreas.wilm@gmail.com>
# License: WTFPL http://www.wtfpl.net/

# Determine number of reads contained in fastq file (given as
# argument; possibly zipped). if a corresponding (newer) fastqc file
# can be found use that instead.
#
# There are many ways to count the number of reads. In most cases you
# can just divide the number of lines by four, however the spec
# does not require the sequence itself or the quality to be on just
# one line. Furthermore zgrep -v '^@' might return hits on
# base-quality lines as well and restricting it further is not
# possible since the fasta-id is optional. The following also assumes
# that a read spans four lines, but it does at least perform some
# checks and supports gzip.
#
# Could check every fourth line in awk (if (NR%4 != 1) {next}) but
# that takes surprisingly longer than getting every fourth line with
# sed first (using non-gnu sed call for compatibility here. The
# following would work in gnu sed as well: -n '1~4p')

for fastq in "$@"; do
    echo -ne "$fastq\t"
    fastqc=${fastq%.gz}
    fastqc=${fastqc%.*}_fastqc/fastqc_data.txt
    if [ -e "$fastqc" ] && [ "$fastqc" -nt $"fastq" ] ; then
        grep '^Total Sequences' "$fastqc" | cut -f 2
    else
        { zcat $fastq 2>/dev/null || cat $fastq; } | \
            #sed -n '1~4p' | \
            sed -n '1,${p;n;n;n;}' | \
            awk '{if (substr($0,0,1)!="@") {exit 1}; nreads+=1;} END {print nreads}'
    fi
done
