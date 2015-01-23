#!/usr/bin/env python
"""Print variants in at least two of given vcf files. Output is not sorted
"""

import sys
import vcf
from collections import Counter



def key_to_var(key):
    """decodes key generated with var_key() and returns it as a generic vcf record"""
    (chrom, pos, ref, alt) = key.split("-")
    return '\t'.join([chrom, pos, ".", ref, alt, ".", "PASS", "."])


def var_key(var):
    """generate a key for a (pyvcf) variant"""


    if len(var.ALT)>1:# FIXME
        sys.stderr.write("Multi allelic entry found. Keeping ony first of %s\n" % var)
    # .ALT is not as string
    alt = str(var.ALT[0])
    return "%s-%s-%s-%s" % (
        var.CHROM, var.POS, var.REF, alt)


if __name__ == "__main__":
    VCFS = sys.argv[1:]
    assert len(VCFS)>1, ("Expected at least two vcf files as argument")
    
    all_var_keys = []
    for f in VCFS:
        vcfr = vcf.Reader(filename=f)
        for v in vcfr:
            if not v.FILTER:
                all_var_keys.append(var_key(v))

    # plain vcf output. creating new record for vcf.Writer would be more difficult
    print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    ctr = Counter(all_var_keys)
    for k in [k for (k, v) in ctr.iteritems() if v > 1]:
        print key_to_var(k)
    
