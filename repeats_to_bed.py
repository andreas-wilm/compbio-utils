#!/usr/bin/env python

# Reads repeat masked (repeats converted to Ns) fasta file and prints
# coordinates of repeat regions in bed format

import os
import sys
from itertools import groupby
from Bio import SeqIO

assert len(sys.argv)==2, (
    "ERROR: Need fasta file as only argument (with repeats marked as N)\n")
fasta=sys.argv[1]
#fasta = 'Homo_sapiens_assembly19.fasta'
assert os.path.exists(fasta)


fh = open(fasta)
seqrec_it = SeqIO.parse(fh, 'fasta')
# debugging
#seqrec_it = ["NACNAGCAGCGAGCAGGCAGNNNNNAGCAGCGANNNNAGCAGCNANNGCANN"]

for seqrec in seqrec_it:
    block_start = 0
    for (is_n, grp) in groupby(seqrec.seq, lambda x: x=='N'):
        # only keep length of stretch (if not in debugging mode)
        l = len(list(grp))
        
        # or instead debug: create string
        #s = ''.join(grp)
        #l = len(s)
            
        if is_n:
            # print bed
            print "%s\t%d\t%d" % (seqrec.id, block_start, block_start+l)
        #else:
        #    # debug
        #    #print seqrec.id, block_start, len(s), ":", s[:5], "...", s[-5:], "==", seqrec.seq[block_start:block_start+5]

        block_start += l
 
fh.close()

























