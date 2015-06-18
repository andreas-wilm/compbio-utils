#!/usr/bin/env python
"""Converts variant positions in given annotated vcf (SNPeff) into
 relative CDS positions
"""

import sys
from math import ceil

import vcf

def main(ann_vcf, scale=1000):
    """main function: takes annotated vcf file as input and spits out 
    relative variant positions
    """

    print "#CHROM\tPOS\tRELPOS(MAX={})".format(scale)

    if ann_vcf == "-":
        vcfr = vcf.VCFReader(sys.stdin)
    else:
        vcfr = vcf.VCFReader(filename=ann_vcf)

    try:
        cds_idx = dict(vcfr.infos)['ANN'].desc.split('|').index(
            ' CDS.pos / CDS.length ')
    except:
        sys.stderr.write(
            "CDS annotation in ANN INFO field of {} not found".format(ann_vcf))
        raise
    for var in vcfr:
        if not var.INFO.has_key('ANN'):
            sys.stderr.write(
                "WARNING: {}:{} has no anntation key in INFO field. Skipping...\n".format(
                    var.CHROM, var.POS, ann_vcf))
            continue
        for ann in var.INFO['ANN']:
            cds_info = ann.split('|')[cds_idx]
            if not cds_info:
                if not 'intergenic_region' in ann:
                    sys.stderr.write(
                        "No CDS info found for var {}:{} but not intergenic\n".format(
                            var.CHROM, var.POS))
            else:
                pos, length = [int(x) for x in cds_info.split('/')]
                fake_pos = int(ceil(pos / float(length) * scale))
                print "{}\t{}\t{}".format(var.CHROM, var.POS, fake_pos)
                # only take first matching annotation content into 
                # account for each var:
                break
    
if __name__ == "__main__":
    main(sys.argv[1])
