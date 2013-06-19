#!/usr/bin/env python
"""Prune certain columns from alignment
"""


#--- standard library imports
#
import sys
import logging
# optparse deprecated from Python 2.7 on
from optparse import OptionParser, SUPPRESS_HELP
from collections import Counter


#--- third-party imports
#
import Bio
from Bio import AlignIO

#--- project specific imports
#
import bioutils

                                                        
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
    creates an OptionParser instance
    """

    # http://docs.python.org/library/optparse.html
    usage = "%prog: " + __doc__ + "\n" \
            "usage: %prog [options]"
    parser = OptionParser(usage=usage)

    parser.add_option("", "--verbose",
                      dest="verbose",
                      action="store_true",
                      help=SUPPRESS_HELP) #"be verbose")
    parser.add_option("", "--debug",
                      dest="debug",
                      action="store_true", 
                      help=SUPPRESS_HELP) #"debugging")
    parser.add_option("", "--all-gap",
                      action="store_true", 
                      dest="all_gap",
                      help="Prune columns if all residues are gaps")
    parser.add_option("", "--any-gap",
                      action="store_true", 
                      dest="any_gap",
                      help="Prune columns with at least one gap")
    parser.add_option("", "--identical",
                      action="store_true", 
                      dest="identical",
                      help="Prune columns if all residues are identical")
    parser.add_option("-i", "--in",
                      dest="aln_in",
                      help="Input alignment ('-' for stdin)")
    parser.add_option("-f", "--infmt",
                      dest="informat",
                      help="Input format (must be supported by Biopython)")
    return parser


def prune_aln(aln, what, fh_out=sys.stdout):
    """Prune what columns from alignment and print result
    """

    assert what in ['any_gap', 'all_gap', 'identical']
    
    keep_cols = []
    for i in xrange(aln.get_alignment_length()):
        # deprecated: col = aln.get_column(i)
        col_nucs = [sr.seq[i].upper() for sr in aln]
        counter = Counter(col_nucs)

        if what == 'any_gap':
            if any([bioutils.isgap(c) for c in counter.keys()]):
                continue
        if what == 'all_gap':
            if all([bioutils.isgap(c) for c in counter.keys()]):
                continue
        if what == 'identical':
            if len(set(counter.keys())) == 1:
                continue

        keep_cols.append(i)

    # FIXME add support for proper alignment output, not just
    # concatenated fasta
    LOG.info("Keeping the following columns: %s" % (
        ', '.join([str(x+1) for x in keep_cols])))
    for s in aln:
        fh_out.write(">%s\n" % s.id)
        fh_out.write('%s\n' % ''.join([s.seq[i] for i in keep_cols]))


        
def main():
    """
    The main function
    """

    parser = cmdline_parser()
    (opts, args) = parser.parse_args()

    if opts.verbose:
        LOG.setLevel(logging.INFO)
    if opts.debug:
        LOG.setLevel(logging.DEBUG)
        
    if not opts.aln_in:
        parser.error("Missing input alignment argument")
        sys.exit(1)

    what = None
    if opts.any_gap:
        assert not what, ("Can only do one operation at a time")
        what = 'any_gap'
    if opts.all_gap:
        assert not what, ("Can only do one operation at a time")
        what = 'all_gap'
    if opts.identical:
        assert not what, ("Can only do one operation at a time")
        what = 'identical'
    if not what:
        parser.error("No operation selected")
        sys.exit(1)
        
    if opts.aln_in == "-":
        fh_in = sys.stdin
    else:
        fh_in = open(opts.aln_in, "rU")

    fmt = opts.informat
    if not fmt:
        fmt = bioutils.guess_seqformat(opts.aln_in)

    aln = AlignIO.read(fh_in, fmt)
    if fh_in != sys.stdin:
        fh_in.close()


    prune_aln(aln, what, sys.stdout)

    

if __name__ == "__main__":
    
    main()

    if sys.version_info < (2 , 7) or sys.version_info > (2 , 8):
        LOG.info("only tested Python 2.7 so far")

    biopython_version = tuple(
        [int(x) for x in Bio.__version__.split('.')])
    if biopython_version < (1 , 55) or biopython_version > (1 , 57):
        LOG.info("using untested version of Biopython")
    
    LOG.info("Successful exit")
