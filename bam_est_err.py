#!/usr/bin/env python
"""Naive error rate estimation from mapping which is assumed to be clonal"""

import sys
import pysam# pre 0.8

BAM_CMATCH = 0
BAM_CINS = 1
BAM_CDEL = 2

sam = pysam.Samfile(sys.argv[1])
opcount = dict()
for read in sam:
    if read.is_unmapped or read.is_duplicate or read.is_qcfail:
        continue
    for (op, ln) in read.cigar:
        if op in [0, 1, 2]:
            opcount[op] = opcount.get(op, 0) + ln

print "#(mis)matches=%d" % opcount[BAM_CMATCH]
print "%%insertions=%.2f" % (opcount[BAM_CINS]/float(opcount[BAM_CMATCH])*100)
print "%%deletions=%.2f" % (opcount[BAM_CDEL]/float(opcount[BAM_CMATCH])*100)
