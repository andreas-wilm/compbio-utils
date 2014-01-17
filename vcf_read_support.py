#!/usr/bin/env python
"""Extract reads supporting variants (or reference) listed in vcf file

Output SAM with extra tag VV:Z:chr:pos-pos:ref>alt
where chr, pos, ref and alt correspond to the variant of question
"""

__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2014 Genome Institute of Singapore"
__license__ = "GPL2"


#--- standard library imports
#
import sys
import logging
import os
import argparse
import gzip
from collections import namedtuple

#--- third-party imports
#
import pysam


#--- project specific imports
#
# /


#global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')


Variant = namedtuple('Variant', 
    ['chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info'])
# all fields are strings with the exception of:
#   pos: int (-1 based)
#   qual: an int if not missing, otherwise "."
#   info: dict


SKIP_FLAGS = [0x4, 0x100, 0x200, 0x400]


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
    parser.add_argument("-b", "--bam",
                        required=True,
                        help="Input BAM file matching vcf")
    parser.add_argument("-i", "--vcf",
                        required=True,
                        help="Input VCF file containing variants to analyze")
    default = 0
    parser.add_argument("--mq-filter",
                        dest="min_mq",
                        type=int,
                        default=default,
                        help="Ignore reads with mapping quality below this value (default=%d)" % default)
    default = 5
    parser.add_argument("--bq-filter",
                        dest="min_bq",
                        type=int,
                        default=default,
                        help="Ignore reads with bases below this value (default=%d)" % default)
    parser.add_argument("-a", "--use-orphan",
                        action="store_true",
                        help="Don't ignore orphan-reads / anomalous read-pairs")

    parser.add_argument("-r", "--ref-only",
                        action="store_true",
                        help="Print reads supporting reference"
                        " instead of variant base")

    return parser



def simple_vcf_reader(fh):
    """yields Variant (chromosome and position only) for variants
    listed in vcf file
    """

    for line in fh:
        if line.startswith('#'):
            continue
        ls = line.rstrip().split('\t')
        # 8 fixed fields per record
        assert len(ls)>=8, (
            "Number of retrieved fields in vcf file too small")
        # ignoring the rest
        (chrom, pos, id, ref, alt, qual, filter, info) = ls[:8]
        pos = int(pos)-1
        try:
            qual = int(qual)
        except:
            qual = "."
        info = dict([field.split('=') for field in info.split(';')])
        yield Variant(chrom, pos, id, ref, alt, qual, filter, info)
        
        
def main():
    """The main function
    """
    
    parser = cmdline_parser()
    args = parser.parse_args()
    
    if args.verbose:
        LOG.setLevel(logging.INFO)
    if args.debug:
        LOG.setLevel(logging.DEBUG)
        import pdb
        from IPython.core import ultratb
        sys.excepthook = ultratb.FormattedTB(mode='Verbose',
                                             color_scheme='Linux', call_pdb=1)        
    
    assert os.path.exists(args.bam), (
        "BAM file %s does not exist" % args.bam)
    sam_in_fh = pysam.Samfile(args.bam)

    sam_out_fh = pysam.Samfile("-", "w", template=sam_in_fh)

    # setup vcf_reader
    # 
    if args.vcf == '-':
        vcf_reader = simple_vcf_reader(sys.stdin)
    else:
        if args.vcf[-3:] == '.gz':
            vcf_reader = simple_vcf_reader(gzip.open(args.vcf))
        else:
            vcf_reader = simple_vcf_reader(open(args.vcf))
            
    variants = [r for r in vcf_reader]
    LOG.info("Loaded %d variants from %s" % (len(variants), args.vcf))
    
    for var in variants:
        if var.info.has_key('INDEL'):
            LOG.warn("Skipping unsupported indel variant at %s:%d" % (
                var.chrom, var.pos+1))
            continue
        if len(var.ref)>1 or len(var.alt)>1:
            LOG.warn("Skipping ref/alt variant with more than"
                     " 1 base at %s:%d" % (var.chrom, var.pos+1))
            continue
            
        reads = list(sam_in_fh.fetch(reference=var.chrom,
                                 start=var.pos, end=var.pos+1))
        LOG.info("%s %d: %d (unfiltered) reads covering position" % (
           var.chrom, var.pos+1, len(reads)))

        for r in reads:

            # FIXME combine
            for f in SKIP_FLAGS:
                if r.flag & f:
                    continue
                
            orphan = (r.flag & 0x1) and not (r.flag & 0x2)
            if orphan and not args.use_orphan:
                continue

            if r.mapq < args.min_mq:
                continue
        
            vpos_on_read = [vpos_on_read 
                            for (vpos_on_read, vpos_on_ref) in r.aligned_pairs 
                            if vpos_on_ref==var.pos]
            assert len(vpos_on_read)==1
            vpos_on_read = vpos_on_read[0]
            if vpos_on_read == None:# skip deletions
                continue

            b = r.query[vpos_on_read]
            bq = ord(r.qqual[vpos_on_read])-33

            if bq < args.min_bq:
                continue

            has_ref = False
            has_var = False
            if b.upper() == var.ref[0].upper():
                has_ref = True
            elif b.upper() == var.alt[0].upper():
                has_var = True
            else:
                # ignore non ref non var
                continue

            # only way I found to add tags. inspired by
            # http://www.ngcrawford.com/2012/04/17/python-adding-read-group-rg-tags-to-bam-or-sam-files/
            var_tag = ('VV', '%s:%d-%d:%s>%s' % (
                var.chrom, var.pos+1, var.pos+1, var.ref, var.alt))
            new_tags = r.tags
            new_tags.append(var_tag)
            r.tags = new_tags
             
            if has_ref and args.ref_only:
                sam_out_fh.write(r)
            elif has_var and not args.ref_only:
                sam_out_fh.write(r)

    sam_in_fh.close()
    # FIXME close sam out if not stdout
    
    
if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
