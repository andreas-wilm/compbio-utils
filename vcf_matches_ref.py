#!/usr/bin/env python3
"""check whether contig in vcf all match fasta entries
"""

# standard library imports
import sys
import os

# third-party imports
import vcf
from pyfaidx import Fasta

# project specific import
# /


def main(vcffile, fasta):
    """Main function
    """
    vcfreader = vcf.VCFReader(filename=vcffile)
    refseqs = Fasta(fasta)

    vcf_contig_to_length_map = dict(
        (contig.id, contig.length) for contig in vcfreader.contigs.values())

    #print("DEBUG: vcf_contig_to_length_map={}".format(vcf_contig_to_length_map))
    
    if len(vcf_contig_to_length_map)==0:
        sys.stderr.write("WARN: no contigs found in vcf\n")
        sys.exit(1)
    num_sq_checked = len(vcf_contig_to_length_map)
    
    num_missing_sq = num_len_mismatch = 0
    for sq in refseqs:
        sq_name, sq_len = sq.name, len(sq)
        if sq_name not in vcf_contig_to_length_map:
            sys.stderr.write("WARN: Fasta sequence {} not found in vcf\n".format(sq_name))
            num_missing_sq += 1
        else:
            if sq_len != vcf_contig_to_length_map[sq.name]:
                sys.stderr.write("WARN: Length mismatch for sequence {}\n".format(sq_name))
                num_len_mismatch += 1

    print("Number of SQs checked (src:vcf): {}".format(num_sq_checked))
    print("Missing SQs (vcf->fasta): {}".format(num_missing_sq))
    print("Lengh mismatches: {}".format(num_len_mismatch))

    # FIXME not checking the other way around




if __name__ == "__main__":

    try:
        vcffile = sys.argv[1]
        fasta = sys.argv[2]
    except IndexError:
        myname = os.path.basename(sys.argv[0])
        sys.stderr.write("Usage: {} {} {}\n".format(myname))
        sys.exit(1)
    main(vcffile, fasta)
