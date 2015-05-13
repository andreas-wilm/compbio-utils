#!/usr/bin/env python
"""Infer indel AF from Strelka VCF

See also 
- https://sites.google.com/site/strelkasomaticvariantcaller/home/somatic-variant-output
- https://groups.google.com/forum/#!msg/strelka-discuss/NAxYqCys4Gg/5p7liDwrXPgJ
- https://groups.google.com/forum/#!topic/strelka-discuss/g_Muy5wVjbY
"""


# --- standard library imports
import os, sys


# --- third party import
import vcf


def strelka_indel_af(vcf_file):
    """Print basic info for each indel variant in strelka vcf and adds
    indel AF for each sample
    """

    if vcf_file == "-":
        vcfreader = vcf.VCFReader(sys.stdin)
    else:
        assert os.path.exists(vcf_file)
        vcfreader = vcf.VCFReader(filename=vcf_file)
        
    print "CHROM\tPOS\t{}".format('\t'.join(vcfreader.samples))
    for var in vcfreader:
        assert var.is_indel
        # print minimal variant info
        print "{}\t{}".format(var.CHROM, var.POS),
        for s in range(len(var.samples)):
            tar = [int(x) for x in var.samples[s].data.TAR]
            tir = [int(x) for x in var.samples[s].data.TIR]
            #print "tar", tar, " tir", tir
            tar = sum(tar)
            tir = sum(tir)
            print "\t{}".format(tir/float(tir+tar)),
        print


if __name__ == "__main__":
    assert len(sys.argv)==2, ("Need strelka vcf as input")
    vcf_file = sys.argv[1]
    strelka_indel_af(vcf_file)
