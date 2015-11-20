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

    print "#CHROM\tPOS\tRELPOS(MAX={})\tIMPACT\tGENE".format(scale)

    if ann_vcf == "-":
        vcfr = vcf.VCFReader(sys.stdin)
    else:
        vcfr = vcf.VCFReader(filename=ann_vcf)


    idx = dict()
    for (key, pattern) in [('cds', ' CDS.pos / CDS.length '),
                            ('impact', ' Annotation_Impact '),
                            ('gene', ' Gene_ID ')]:
        try:
            idx[key] = dict(vcfr.infos)['ANN'].desc.split('|').index(pattern)
            #sys.stderr.write("DEBUG: idx[{}]={}\n".format(key, idx[key]))
        except:
            sys.stderr.write(
                "{} annotation in ANN INFO field of {} not found".format(key, ann_vcf))
            raise
                   
    for var in vcfr:
        if not var.INFO.has_key('ANN'):
            sys.stderr.write(
                "WARNING: {}:{} has no annotation key in INFO field. Skipping...\n".format(
                    var.CHROM, var.POS, ann_vcf))
            continue
        #if len(var.INFO['ANN'])>1: sys.stderr.write("DEBUG: >1 ANN for {}".format(var))
        # printed once per annotation
        for ann in var.INFO['ANN']:
            info = dict()
            for k in idx.keys():
                info[k] = ann.split('|')[idx[k]]
                #sys.stderr.write("DEBUG: info[{}]={}\n".format(k, info[k]))
            if not info['cds']:
                if not 'intergenic_region' in ann:
                    sys.stderr.write(
                        "No CDS info found for var {}:{} but not intergenic\n".format(
                            var.CHROM, var.POS))
            else:
                pos, length = [int(x) for x in info['cds'].split('/')]
                fake_pos = int(ceil(pos / float(length) * scale))
                print "{}\t{}\t{}\t{}\t{}".format(var.CHROM, var.POS, fake_pos, info['impact'], info['gene'])
    
if __name__ == "__main__":
    main(sys.argv[1])
