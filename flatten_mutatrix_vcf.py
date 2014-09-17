#!/usr/bin/env python
"""Flatten mutatrix vcf into separate entries, i.e. resolving multi-allelic events
"""

__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2014 Genome Institute of Singapore"
__license__ = "WTFPL http://www.wtfpl.net/"


# --- standard library imports
#
import sys
import os

#--- third-party imports
#
import vcf


def print_snp(stream, chrom, pos, ref, alt, qual, info):
    """simplified printing of variant
    """
    stream.write("%s\t%s\t.\t%s\t%s\t%s\t.\t%s\n" % (
            chrom, pos, ref, alt, qual, info))


def print_vcf_header(stream):
    """print vcf header
    """
    stream.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")


def main():
    """main function (shutup pylint)
    """

    assert len(sys.argv)==2
    f = sys.argv[1]
    assert os.path.exists(f)

    print_vcf_header(sys.stdout)

    vcfreader = vcf.VCFReader(filename=f)
    for v in vcfreader:
        assert len(v.ALT) == len(v.INFO['TYPE'])

        for i in range(len(v.ALT)):
            t = v.INFO['TYPE'][i]
            a = str(v.ALT[i])
            if t == 'snp':
                print_snp(sys.stdout, v.CHROM, v.POS, v.REF, a, v.QUAL, "snp")
            elif t == 'mnp':
                assert len(v.REF)>1
                assert len(v.REF)==len(a)
                for i in range(len(v.REF)):
                    print_snp(sys.stdout, v.CHROM, v.POS+i, 
                              v.REF[i], a[i], v.QUAL, "mnp")
            else:
                # FIXME handle indels
                pass
    

if __name__ == "__main__":
    main()
            
            
