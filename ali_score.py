#!/usr/bin/env python
"""Score an alignment

Implement gap handling and more scores:
http://www.biostars.org/post/show/3856/entropy-from-a-multiple-sequence-alignment-with-gaps/

William Valdar : Scoring Residue Conservation.
http://onlinelibrary.wiley.com/doi/10.1002/prot.10146/abstract;jsessionid=18D6E98A259624E3D6616386C0EC32C5.d03t02
"""


#--- standard library imports
#
import os
import sys
from math import log
import logging
# optparse deprecated from Python 2.7 on
from optparse import OptionParser
from collections import Counter

#--- third-party imports
#
from Bio import AlignIO

#--- project specific imports
#
import bioutils

                                                        
__author__ = "Andreas Wilm"
__version__ = "0.1"
__email__ = "andreas.wilm@gmail.com"
__license__ = "The MIT License (MIT)"


#global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')




def shannon_entropy(l, b=2):
    """Return the Shannon entropy of random variable with probability
    vector l.

    Adopted from
    http://www.lysator.liu.se/~jc/mthesis/A_Source_code.html#functiondef:entropies.py:shannon_entropy
    """
    return sum([-p*log(p, b) for p in l if p > 0])

    
def cmdline_parser():
    """
    creates an OptionParser instance
    """

    # http://docs.python.org/library/optparse.html
    usage = "%prog: " + __doc__ + "\n" \
            "usage: %prog [options]"
    parser = OptionParser(usage=usage)

    parser.add_option("", "--verbose",
                      action="store_true", dest="verbose",
                      help="be verbose")
    parser.add_option("", "--debug",
                      action="store_true", dest="debug",
                      help="debugging")

    parser.add_option("-i", "--aln",
                      dest="aln_in",
                      default="-",
                      help="Input alignment or '-' for stdin (default)")
    return parser


def seqid(colctr):
    """FIXME:add-doc
    """

    count = 0
    # sort colctr counter dictionary in descending order and ignore gaps
    for (res, count) in sorted(colctr.items(), key=lambda x: x[1], reverse=True):
        if res in "-~.":
            continue
        break
    # res is now the most frequent residues, with count counts
    #import pdb; pdb.set_trace()
    
    return count/float(sum(colctr.values()))
    

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
        parser.error("Missing input alignment argument\n")
        sys.exit(1)
    if len(args):
        parser.error("Unrecognized arguments found: %s" % args)

    if False:
        char_set = "ACGTU"
    else:
        char_set = "ACDEFGHIKLMNPQRSTVWY"
        #x = any
        #z = Gln or Glu
        #b = Asp or Asn
    LOG.warn("hardcoded charset %s" % char_set)
    
    if opts.aln_in != "-" and not os.path.exists(opts.aln_in):
        LOG.fatal("Input alignment %s does not exist.\n" % opts.aln_in)
        sys.exit(1)

    if opts.aln_in == "-":
        fh = sys.stdin
        fmt = 'fasta'
    else:
        fmt = bioutils.guess_seqformat(opts.aln_in)
        fh = open(opts.aln_in, "rU")
                
    entropy_per_col = []    
    seqid_per_col = []    
    aln = AlignIO.read(fh, fmt)
    for i in xrange(aln.get_alignment_length()):
        col = aln.get_column(i)
        counter = Counter(col.upper())

        vec = []
        # this will ignore invalid chars incl. ambiguities        
        denom = sum([counter[r] for r in char_set])
        if denom == 0:
            LOG.fatal("denom = 0, means no valid chars in col %d?" % (i+1))
            #import pdb; pdb.set_trace()
            raise ValueError
        for res in char_set:
            vec.append(counter[res]/float(denom))
        LOG.debug("vec=%s denom=%s counter=%s" % (vec, denom, counter))
        entropy_per_col.append(shannon_entropy(vec))

        seqid_per_col.append(seqid(counter))
        
        print "%d %.6f %.6f %s" % (
            i+1, seqid_per_col[i], 
            entropy_per_col[i], 
            ' '.join(["%s:%d" % (k,v) 
                      for (k,v) in sorted(counter.iteritems())]))

    if fh.name != '-':
        fh.close()

    #for i in range(len(entropy_per_col)):
    #print i+1, entropy_per_col[i]

        
if __name__ == "__main__":
    main()


