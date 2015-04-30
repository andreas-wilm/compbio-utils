#!/usr/bin/env python
"""Extract values from vcf file
"""

# --- standard library imports
import os, sys
import argparse

# --- third party import
import vcf

def main():
    """main function
    """
        
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("-i", "--vcf",
                        required="True",
                        dest="vcf",
                        help="VCF file to read (- for stdin)")
    parser.add_argument("-a", "--all",
                        dest="ign_filter",
                        help="Use all, not just passed variants")
    choices = ['SNV', 'INDEL']
    parser.add_argument("-t", "--type",
                        dest="type",
                        choices=choices,
                        help="Only work on this type of variants (one of %s)" % (','.join(choices)))
    parser.add_argument("-v", "--value",
                        dest="value",
                        required=True,
                        help="What to print. Must be 'QUAL'"
                        " or any field name from vcf INFO.")
    args = parser.parse_args()

    if args.value == "QUAL":
        def extract_func(var):
            """Extract quality value taken case of missing values"""
            if not var.QUAL and  var.QUAL != 0:
                return "."
            else:
                return var.QUAL
    else:
        def extract_func(var):
            """Extract value from vcf INFO field"""
            val = var.INFO[args.value]
            if isinstance(val, list):
                return ','.join([str(x) for x in var.INFO[args.value]])
            else:
                return val

    if args.vcf == "-":
        vcfreader = vcf.VCFReader(sys.stdin)
    else:
        assert os.path.exists(args.vcf)
        vcfreader = vcf.VCFReader(filename=args.vcf)        
        
    for var in vcfreader:
        if not args.ign_filter and var.FILTER:
            #print "Skipping filtered %s" % str(var)
            continue
        if args.type == 'INDEL' and not var.is_indel:
            #print "Skipping non-INDEL %s" % var
            continue
        if args.type == 'SNV' and not var.is_snp:
            #print "Skipping non-SNV %s " % var
            continue
        print extract_func(var)

if __name__ == "__main__":
    main()
        
