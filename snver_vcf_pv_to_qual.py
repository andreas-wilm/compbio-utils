#!/bin/env python
"""add quality field to snver vcf by converting pvalue in info field
"""

import sys
from lofreq import utils

# - pv in snver raw is same as in snver filtered
# - no additional filters are applied at that stage AFAIK
# - cleaner output can be achiever by removing anything with Q<1 and not PASSED

vcf = sys.argv[1]
with open(vcf) as fh:
    for line in fh:
        if line.startswith('#'):
            print line,
        else:
            pv = float(line[line.index('PV=')+3:].split()[0])
            q = utils.prob_to_phredqual(pv)
            line = line.replace('.\tPASS', '%d\tPASS' % q)
            print line,

