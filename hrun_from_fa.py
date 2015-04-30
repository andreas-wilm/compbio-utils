#!/usr/bin/env python
"""Extract homopolymer runs from fasta file and save as bed
"""

import os
import sys
from itertools import groupby

import pysam


def hrun_from_fa(fasta, bed):
    if bed == "-":
        bedfh = sys.stdout
    else:
        assert not os.path.exists(bed)
        bedfh = open(bed, 'w')
    fafh = pysam.Fastafile(fasta)

    for ref in fafh.references:
        pos = 0
        # no case conversion. always requires copy of seq in mem and might
        # not be wanted in th first place.
        for (k, g) in groupby(fafh.fetch(ref)):
            hrun = ''.join(g)
            if len(hrun)>=2:
                bedfh.write("{}\t{:d}\t{:d}\tHRUN-{:s}\t{:d}\n".format(
                        ref, pos, pos+len(hrun), k, len(hrun)))
            pos += len(hrun)
    if bedfh != sys.stdout:
        bedfh.close()


if __name__ == "__main__":
    try:
        fasta = sys.argv[1]
        bed = sys.argv[2]
    except IndexError:
        sys.stderr.write("FATAL: fasta and bed file argument are required\n")
        sys.exit(1)
    
    hrun_from_fa(fasta, bed)
    
