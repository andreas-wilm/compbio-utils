#!/usr/bin/env python
"""Get freqs for nucleotide pairs (present per read) at given
positions."""


#--- standard library imports
#
import sys
import os
import logging
import argparse

#--- third-party imports
#
import pysam


#--- project specific imports
#
# /


__author__ = "Andreas Wilm"
__version__ = "0.1"
__email__ = "andreas.wilm@gmail.com"
__license__ = "The MIT License (MIT)"


# global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN, 
                    format='%(levelname)s [%(asctime)s]: %(message)s')



def cmdline_parser():
    """
    creates an argparse instance
    """

    # http://docs.python.org/library/optparse.html
    parser = argparse.ArgumentParser(
        description=__doc__)

    parser.add_argument("--verbose",
                      dest="verbose",
                      action="store_true",
                      help=argparse.SUPPRESS) #"be verbose")
    parser.add_argument("--debug",
                      dest="debug",
                      action="store_true", 
                      help=argparse.SUPPRESS) #"debugging")
    parser.add_argument("-1", "--pos1",
                      dest="pos1",
                      type=int,
                      required=True,
                      help="First position")
    parser.add_argument("-2", "--pos2",
                      dest="pos2",
                      type=int,
                      required=True,
                      help="Second position")
    parser.add_argument("-r", "--ref",
                      dest="ref",
                      required=True,
                      help="Mapping/reference sequence/chromosome name")
    parser.add_argument("-b", "--bam",
                      dest="bam",
                      required=True,
                      help="Mapping input file (BAM)")
    parser.add_argument("--fasta",
                        dest="fasta",
                        help="Print corresponding reference nucleotides"
                        " region from this fasta file as well")
    return parser



def main():
    """
    The main function
    """

    parser = cmdline_parser()
    args = parser.parse_args()
    
    if args.verbose:
        LOG.setLevel(logging.INFO)
    if args.debug:
        LOG.setLevel(logging.DEBUG)
        
    if args.pos1 >= args.pos2:
        LOG.fatal("First position must be smaller than second")
        parser.print_usage(sys.stderr)
        
    # file check
    if not os.path.exists(args.bam):
        LOG.fatal("file '%s' does not exist.\n" % args.bam)
        sys.exit(1)
    if args.fasta and not os.path.exists(args.fasta):
        LOG.fatal("file '%s' does not exist.\n" % args.fasta)
        sys.exit(1)
        
    pos_pair = (args.pos1-1, args.pos2-1)
    
    sam = pysam.Samfile(args.bam, "rb")

    if args.fasta:
        fastafile = pysam.Fastafile(args.fasta)
        refregion = fastafile.fetch(args.ref, pos_pair[0], pos_pair[1]+1)
        print "Ref. region: %s" % refregion
    
    # initialize counts to valid nucleotide combinations even though
    # not strictly necessary, since we deal with non existing values
    # later.
    #
    VALID_NUCS = 'ACGTN'
    counts = dict()
    for n1 in VALID_NUCS:
        # counting twice. booooooh
        for n2 in VALID_NUCS:
            counts["%c%c" % (n1, n2)] = 0
            
    num_dups = 0
    num_anomalous = 0
    for alnread in sam.fetch(args.ref, pos_pair[0], pos_pair[1]+1):
        if alnread.is_duplicate:
            num_dups += 1
            continue
        assert not alnread.is_unmapped # shouldn't be possible anyway
        if alnread.is_paired and not alnread.is_proper_pair:
            num_anomalous += 1
            continue
        #region = list(sam.fetch("gi|158976983|ref|NC_001474.2|", 10308, 10411))
        #[(i, str(r)) for (i, r) in enumerate(region) if r.cigar != [(0, 51)]]

        # create a map of ref position (string as key) and the corresponding query
        # (clipped read) nucleotides (value)
        pos_nt_map = dict([(str(rpos), alnread.query[qpos])
                           for (qpos, rpos) in alnread.aligned_pairs 
                           if qpos and rpos])

        try:
            nt5 = pos_nt_map[str(pos_pair[0])]
            nt3 = pos_nt_map[str(pos_pair[1])]
        except KeyError:
            #print "read not fully overlapping both positions (want %d-%d, have %d-%d, cigar %s): %s" % (
            #    pos_pair[0]+1, pos_pair[1]+1, 
            #    alnread.aend-alnread.alen, alnread.aend, 
            #    alnread.cigar, alnread.positions)
            continue
        key = "%c%c" % (nt5, nt3)
        counts[key] = counts.get(key, 0) + 1
        #if key == 'AT':
        #    import pdb; pdb.set_trace()
        
        # we've initialized counts but be paranoid about existance of key anyway
    
    if num_dups:
        print "Ignored %d reads flagged as duplicates" % num_dups
    if num_anomalous:
        print "Ignored %d reads flagged as paired but not in proper pair" % num_anomalous
    counts_sum = sum(counts.values())
    print "%d (non-dup)reads overlapped both given positions %d and %d"  % (
        counts_sum, pos_pair[0]+1, pos_pair[1]+1)
    if counts_sum == 0:
        sys.exit(0)
    for k in sorted(counts.keys()):
        if counts[k]:
            print "%s %d %.4f" % (k, counts[k], counts[k]/float(counts_sum))
            
    
    
if __name__ == "__main__":
    main()
    LOG.info("Successful exit")
