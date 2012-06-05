#!/bin/bash

# check if read order in two given paired end reads is the same. done
# via sampling each 10000st read.

fastq_order_md5() {
    zgrep '^@.*/[0-9]$' $1 | \
        awk '{if (NR%10000==0) {sub(/[0-9]$/, ""); print}}' | \
        md5sum | cut -f 1 -d ' '
}


fastq1=$1
fastq2=$2
if [ -z $fastq1 ] || [ -z $fastq2 ]; then
    echo "FATAL: need exactly two arguments: $(basename $0) s1.fastq[.gz] s2.fastq[.gz]" 1>&2
    exit 1
fi
if [ ! -e $fastq1 ] || [ ! -e $fastq2 ]; then
    echo "FATAL: One of the two given files does not exist" 1>&2
    exit 1
fi

md5_s1=$(fastq_order_md5 $fastq1)
md5_s2=$(fastq_order_md5 $fastq2)

if [ $md5_s1 == $md5_s2 ]; then
    echo "OK"
    exit 0
else
    echo "FAILED"
    exit 1
fi

