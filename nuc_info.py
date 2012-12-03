#!/usr/bin/env python
"""FIXME
"""

import sys
from collections import Counter

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_rna

for line in sys.stdin:
    line = line.rstrip()
    if len(line)==0 or line.startswith("#"):
        continue
    
    line = line.upper()

    counts = Counter(line)
    n_u = counts.get('U', 0) + counts.get('u', 0)
    n_t = counts.get('T', 0) + counts.get('t', 0)
    if n_u > n_t:
        seq = Seq(line, generic_rna)
        seqtype = "RNA"
    else:
        seq = Seq(line, generic_dna)
        seqtype = "DNA"

    print "Parsed input as: %s" % seq
    print "Seq. type: %s" % seqtype
    print "Length: %d" % len(seq)
    print "Composition: %s" % (' '.join(
            ["%c:%.2f%% " % (nuc, cnt/float(len(seq))) 
             for (nuc, cnt) in counts.iteritems()]))
    print "Complement: %s" % seq.complement()
    print "Reverse complement: %s" % seq.reverse_complement()
    print "Translated: %s" % seq.translate()
    # See http://biopython.org/wiki/Seq
