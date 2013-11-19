#!/bin/bash
# Simple script (no sanity checks are done) to concat fasta sequences,
# so that seqs are on one line (like vienna format). Takes fasta input
# from pipe or file. Main use is quick-and-dirty downstream handling,
# e.g. sequence extraction with "grep -A1 name" or extracting
# subsequences based on position
# awk '{if (/^>/) {print} else {print substr($0, 666, 6)}}' 
# etc etc
# 

# Test correctness:
# cat test.fa | fa_to_vie.sh | sreformat fasta - | md5sum 
# fa_to_vie.sh test.fa | sreformat fasta - | md5sum 
# sreformat fasta test.fa | md5sum 

awk '{
    if (/^>/) {
        if (NR!=1) {
            printf "\n";
        }
        print $0;
    } else {
        printf $0;
    }
}
END {
    printf "\n";
}' $1
