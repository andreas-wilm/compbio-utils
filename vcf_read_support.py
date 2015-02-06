#!/usr/bin/env python
"""Extract reads overlapping variant positions and tag them according
to whether they support a variant or the reference.

Output SAM with extra tag key:Z:chr:pos:ref>alt where chr, pos, ref
and alt correspond to the variant of question. For reads supporting
variants key is 'vv', for those supporting the reference it's 'vr',
otherwise the read will not be written
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
                        help="Input VCF file containing variants to analyze"
                        " (clashes with --var)")
    parser.add_argument("-v", "--var",
                        help="Report reads for this variant only. Format: chr:pos:ref-alt"
                        " (clashes with --vcf)")
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

    return parser



def simple_vcf_reader(fh):
    """yields variants listed in vcf file as 'Variant'
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
        info_d = dict()
        for field in info.split(';'):
            kv = field.split('=')
            # boolean entries get True as value
            if len(kv)==1:
                info_d[kv[0]] = True
            else:
                info_d[kv[0]] = kv[1]
        #try:
        #    info = dict([field.split('=') for field in info.split(';')])
        #except ValueError:
        #    import pdb; pdb.set_trace()
        yield Variant(chrom, pos, id, ref, alt, qual, filter, info_d)

#def simple_vcf_write(var):
#    """write Variant as vcf entry"""
#    c2f = "\t".join([var.chrom, "%d" % (var.pos+1), var.id, var.ref, var.alt, "%s" % (var.qual), var.filter]
#    # FIXME handle True values
#    info = ";".join(["%s=%s" % (k,v) for (k,v) in var.info.items()])
#    return "%s\t%s" % (c2f, info)
                
        
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

    sam_out_fh = pysam.Samfile("-", "wh", template=sam_in_fh)

    # variants
    #
    #
    if args.vcf and args.var:
        LOG.fatal("Please use one: vcf or variant arg, but bot both")
        sys.exit(1)
    if args.vcf: 
        if args.vcf == '-':
            vcf_reader = simple_vcf_reader(sys.stdin)
        else:
            if args.vcf[-3:] == '.gz':
                vcf_reader = simple_vcf_reader(gzip.open(args.vcf))
            else:
                vcf_reader = simple_vcf_reader(open(args.vcf))
        variants = [r for r in vcf_reader]
        LOG.info("Loaded %d variants from %s" % (len(variants), args.vcf))
        
    elif args.var:
        try:
            (chrom, pos, ref_alt) = args.var.split(":")
            (ref, alt) = ref_alt.split('-')
            pos = int(pos)-1
        except:
            LOG.fatal("Couldn't parse variant %s" % args.var)
            sys.exit(1)
        variants = [Variant(chrom, pos, ".", ref, alt, ".", ".", dict())]
        
    else:       
        LOG.critical("Missing vcf or variant argument") 
        sys.exit(1)
        
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
            if vpos_on_read == None:# FIXME no support for deletions
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

            if has_ref:
                var_tag_key = 'vr'
            elif has_var:
                var_tag_key = 'vv'
            else:
                continue# paranoia (already handled above)
            assert var_tag_key not in [t[0] for t in r.tags], (
                "Oops...tag %s already present in read. Refusing to overwrite")
            var_tag = (var_tag_key, '%s:%d:%s>%s' % (
                var.chrom, var.pos+1, var.ref, var.alt))
            new_tags = r.tags
            new_tags.append(var_tag)
            r.tags = new_tags
             
            sam_out_fh.write(r)

    sam_in_fh.close()
    # FIXME close sam out if not stdout

    # FIXME untangle and move to functions
    
    # FIXME add tests:
    #  1:
    #    vcf_read_support.py -b bam -v var | grep -c var
    #    should give same as
    #    vcf_read_support.py -b bam -i vcf -b bam | grep -c var
    # ...
    
if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
