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
        print extract_func(var)

if __name__ == "__main__":
    main()
        
