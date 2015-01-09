#!/usr/bin/env python

import os
import sys
import logging
from collections import namedtuple
import argparse
 
#sys.path.insert(0, '/Users/wilma/GIS/lofreq/lofreq2-gis.git/src/lofreq_python/lofreq_star/')
# PyVCF by default
import vcf

from matplotlib import pyplot as plt
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt


# https://pypi.python.org/pypi/matplotlib-venn
# http://fouryears.eu/2012/10/13/venn-diagrams-in-python/
from matplotlib_venn import venn3, venn3_circles

#global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')


VcfContainer = namedtuple('Vcf', ['name', 'file', 'text', 'vars', 'var_set'])


def cmdline_parser():
    """Returns an argparse instance
    """

    # http://docs.python.org/dev/howto/argparse.html
    parser = argparse.ArgumentParser(description=__doc__)
    
    parser.add_argument("--verbose",
                        action="store_true",
                        help="Be verbose")
    parser.add_argument("--debug",
                        action="store_true",
                        help="Enable debugging")
    parser.add_argument("-o", "--plot",
                        required=True,
                        help="Output plot filename")
    parser.add_argument("-v", "--vcf",
                        required=True,
                        action='append',
                        help="Two to three vcf files in the form vcf:name:[txt]")
    return parser



def var_key(var):
    # pyvcf encodes .alt not as string, ie.
    return "%s-%s-%s-%s" % (
        var.CHROM, var.POS, var.REF, ''.join(["%s" % k for k in var.ALT]))
       



def main():
    """main fcuntion
    """
    
    parser = cmdline_parser()
    args = parser.parse_args()
           
    if args.verbose:
        LOG.setLevel(logging.INFO)
    if args.debug:
        LOG.setLevel(logging.DEBUG)        
    
    if os.path.exists(args.plot):
        LOG.fatal("Cowardly refusing to overwrite already existing file %s" % args.plot)
        sys.exit(1)
        
    variants = []
    
    for arg in args.vcf:
        s = arg.rstrip().split(':')
        assert len(s)==3, ("Excpected argument format vcf:name:[txt] not found")
        (vcffile, name, text) = s
        print "DEBUG vcffile=%s name=%s text=%s" % (vcffile, name, text)
        
        LOG.info("Reading %s" % vcffile)
        # PyVCF can read gzip
        vcfreader = vcf.VCFReader(filename=vcffile)
        vars = [v for v in vcfreader if v.FILTER in [".", "PASS"] or len(v.FILTER)==0]
        LOG.info("Got %d vars" % len(vars))
        var_set = set(dict([(var_key(v), v) for v in vars]))
        variants.append(VcfContainer(name, vcffile, text, vars, var_set))
    
    subsets = [v.var_set for v in variants]
    labels = [v.name + "\n" + v.text for v in variants]
    v = venn3(subsets=subsets, set_labels=labels)
    c = venn3_circles(subsets=subsets, linestyle='solid')
    #plt.show()
    
    plt.savefig(args.plot)

    
if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")    
