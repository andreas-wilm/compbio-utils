#!/usr/bin/env python
"""Prune certain columns from alignment
"""


#--- standard library imports
#
import os
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


GAP_CHARS = ['-', '~', '.']



def isgap(res):
    """Return true if given residue is a gap character
    """
    return (res in GAP_CHARS)



def guess_seqformat(fseq):
    """Guess sequence format from file extension
    """
    default = 'fasta'

    # Table for guessing the alignment format from the file extension. 
    # See http://www.biopython.org/wiki/SeqIO
    #
    # Only define the ones I actually came accors here:
    ext_to_fmt_table = dict(
        aln = 'clustal',
        embl = 'embl',
        fasta = 'fasta',
        fa = 'fasta',
        genbank = 'genbank',
        gb = 'genbank',
        phylip = 'phylip',
        phy = 'phylip',
        ph = 'phylip',
        pir = 'pir',
        stockholm = 'stockholm',
        st = 'stockholm',
        stk = 'stockholm')

    try:
        fext = os.path.splitext(fseq)[1]
        fext = fext[1:].lower()
        fmt =  ext_to_fmt_table[fext]
    except KeyError:
        return default

    return fmt



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
        col = aln.get_column(i)
        counter = Counter(col.upper())

        if what == 'any_gap':
            if any([isgap(c) for c in counter.keys()]):
                continue
        if what == 'all_gap':
            if all([isgap(c) for c in counter.keys()]):
                continue
        if what == 'identical':
            if len(set(counter.keys())) == 1:
                continue

        keep_cols.append(i)

    # FIXME add support for proper alignment output, not just
    # concatenated fasta
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

        
    if opts.aln_in == "-":
        fh_in = sys.stdin
    else:
        fh_in = open(opts.aln_in, "rU")

    fmt = opts.informat
    if not fmt:
        fmt = guess_seqformat(opts.aln_in)

    aln = AlignIO.read(fh_in, fmt)
    if fh_in != sys.stdin:
        fh_in.close()


    prune_aln(aln, what, sys.stdout)

    

if __name__ == "__main__":
    if sys.version_info < (2 , 7):
        sys.stderr.write("WARNING: only tested Python 2.7 so far\n")
    elif sys.version_info > (2 , 8):
        sys.stderr.write("WARNING: only tested Python 2.7 so far\n")

    biopython_version = tuple([int(x) for x in Bio.__version__.split('.')])
    if biopython_version < (1 , 55):
        sys.stderr.write("WARNING: only tested Biopython 1.55 so far\n")
    elif biopython_version > (1 , 55):
        sys.stderr.write("WARNING: only tested Biopython 1.55 so far\n")

    main()
    LOG.info("Successful exit")
