#!/usr/bin/env python
"""Extract values from vcf file
"""

# --- standard library imports
import os
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

    assert os.path.exists(args.vcf), (
        "Non-existing file %s" % args.vcf)
        
    if args.value == "QUAL":
        extract_func = lambda v: v.QUAL
    else:
        def extract_func(var):
            """Extract value from vcf INFO field"""
            val = var.INFO[args.value]
            if isinstance(val, bool):
                return val
            else:
                # iterate in case this was a list. requires bool check above
                return ','.join([str(x) for x in var.INFO[args.value]])

    vcfreader = vcf.VCFReader(filename=args.vcf)
    
    for var in vcfreader:
        print extract_func(var)

if __name__ == "__main__":
    main()
        
