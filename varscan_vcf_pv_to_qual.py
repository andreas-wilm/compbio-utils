#!/bin/env python
"""Add quality to varscan's vcf file
"""

import sys
from lofreq import utils
import gzip

try:
    vcf = sys.argv[1]
    if vcf == "-":
        fh = sys.stdin
    else:
        if vcf.endswith(".gz"):
            fh = gzip.open(vcf)	
        else:
            fh = open(vcf)
except IndexError:
    fh = sys.stdin
    
for line in fh:
    if line.startswith('#'):
        print line,
    else:
        fields = line.split('\t')
        pv_idx = fields[8].split(':').index("PVAL")
        pv = float(fields[9].split(':')[pv_idx])
        q = utils.prob_to_phredqual(pv)
        line = line.replace('.\tPASS', '%d\tPASS' % q)
        
        af_idx = fields[8].split(':').index("FREQ")
        af_str = fields[9].split(':')[af_idx]
        assert af_str[-1] == "%"
        af = float(af_str[:-1])/100.0
        line = line.replace('ADP=', 'AF=%.4f;ADP=' % af)
        
        print line,
if fh != sys.stdin:
    fh.close()
