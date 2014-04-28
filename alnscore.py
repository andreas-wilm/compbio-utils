#!/usr/bin/env python
"""Score an alignment

Implement gap handling and more scores:
http://www.biostars.org/post/show/3856/entropy-from-a-multiple-sequence-alignment-with-gaps/

William Valdar : Scoring Residue Conservation.
http://onlinelibrary.wiley.com/doi/10.1002/prot.10146/abstract;jsessionid=18D6E98A259624E3D6616386C0EC32C5.d03t02

Output is a table with one row per alignment column:
1. pos
2. identity
3. entropy
4. basecounts as base:count tuples
"""


#--- standard library imports
#
import os
import sys
from math import log
import logging
import argparse
from collections import Counter

#--- third-party imports
#
from Bio import AlignIO

#--- project specific imports
#
from bioutils import GAP_CHARS
from bioutils import guess_seqformat
        
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
    """

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("--verbose",
                      action="store_true", dest="verbose",
                      help="be verbose")
    parser.add_argument("--debug",
                      action="store_true", dest="debug",
                      help="debugging")
    parser.add_argument("-i", "--aln",
                      dest="aln_in",
                      help="Input alignment or '-' for stdin")
    parser.add_argument("-m", "--map-to",
                      dest="map_to",
                      help="If given, using unaligned positions of this"
                      " seq instead of aligned pos (skipping pos with"
                      " gaps in this seq)")
    parser.add_argument("--ignore",
                      dest="ign_seqs",
                      nargs="+",
                      help="If given, ignore these sequences for score"
                      " computation")
    return parser


def seqid(colctr):
    """colctr has to be collections.Counter for residues in an
    alignment column. if it contains gaps, they will be considered
    """

    for (res, count) in colctr.most_common():
        return count/float(sum(colctr.values()))
    raise ValueError


#def get_aln_column(aln, col_no):
#    """Replacement for depreacted Alignment.get_column()
#    """
#    assert col_no < aln.get_alignment_length() 
#    # get_column() and get_all_seq() deprecated
#    return ''.join([s.seq[col_no] for s in aln])


def unaln_pos_map(seq):
    """Returns a list of length seq, where each value contains the
    unaligned position (or None if NA)
    """
    map_to_seq_cols = []
    gap_ctr = 0
    for i in range(len(seq)):
        if seq[i] in GAP_CHARS:
            gap_ctr += 1
        map_to_seq_cols.append(i-gap_ctr)
        if map_to_seq_cols[-1] < 0:
            map_to_seq_cols[-1] = None
    return map_to_seq_cols


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
    if not args.aln_in:
        parser.error("Missing input alignment argument\n")
        sys.exit(1)


    #char_set = "ACGTU"
    #char_set = "ACDEFGHIKLMNPQRSTVWY"
    #x = any
    #z = Gln or Glu
    #b = Asp or Asn
    char_set = "ACGTN"
    char_set_ambig = "N"
    
    LOG.warn("using hardcoded charset %s" % char_set)
    # FIXME auto-detection of alphabet)
    
    if args.aln_in != "-" and not os.path.exists(args.aln_in):
        LOG.fatal("Input alignment %s does not exist.\n" % args.aln_in)
        sys.exit(1)

    if args.aln_in == "-":
        fh = sys.stdin
        fmt = 'fasta'
    else:
        fmt = guess_seqformat(args.aln_in)
        fh = open(args.aln_in, "rU")
                
    entropy_per_col = []    
    seqid_per_col = []    
    # note: had one case where this happily read an unaligned file!?
    aln = AlignIO.read(fh, fmt)

    # if requested, get sequence record for the sequence we should
    # positions to
    map_to_seq = None
    if args.map_to:
        map_to_seq = [rec.seq for rec in aln if rec.id == args.map_to]
        if not len(map_to_seq):
            LOG.fatal("Couldn't find a sequence called %s in %s" % (
                args.map_to, fh.name))
            sys.exit(1)
        elif len(map_to_seq)>1:
            LOG.fatal("Find more than one sequence with name %s in %s" % (
                args.map_to, fh.name))
            sys.exit(1)
        map_to_seq = map_to_seq[0]
        map_to_seq_cols = unaln_pos_map(map_to_seq)

    ign_idxs = []
    if args.ign_seqs:
        for s in args.ign_seqs:
            found = False
            for (i, r) in enumerate([r.id for r in aln]):
                if r==s:
                    ign_idxs.append(i)
                    found = True
                    break
            if not found:
                LOG.warn("No match for ignore sequence %s in alignment" % s)
    LOG.debug("ign_idxs = %s" % ign_idxs)

    
    for cidx in xrange(aln.get_alignment_length()):
        col =  list(aln[:, cidx].upper())
        
        # ignore chars as requested from ign_seqs
        for i in ign_idxs:
            del col[i]
        del i

        # replace unknown characters with ambiguity symbol
        unknown_chars = []
        for (i, c) in enumerate(col):
            if c not in char_set and c not in GAP_CHARS:
                unknown_chars.append(c)
                col[i] = char_set_ambig
            elif c in GAP_CHARS:
                col[i] = "-"
                
        unknown_chars = set(unknown_chars)
        if len(unknown_chars):
            LOG.warn("Found unknown characters in col %d (%s) and replaced them with %c" % (
                cidx+1, unknown_chars, char_set_ambig))
            
        counter = Counter(col)
        denom = sum(counter.values())
        if denom == 0:
            LOG.warn("No valid chars in col %d (col=%s)?" % (cidx+1, col))
            #import pdb; pdb.set_trace()
            #raise ValueError
            entropy_per_col.append(-1)
            seqid_per_col.append(-1)
        else:
            vec = []
            # count gaps for entropy 
            for res in list(char_set) + ["-"]:
                vec.append(counter[res]/float(denom))
            LOG.debug("vec=%s denom=%s counter=%s" % (vec, denom, counter))
            entropy_per_col.append(shannon_entropy(vec))
            seqid_per_col.append(seqid(counter))


        # due to the fact that we keep all values (which is actually
        # not necessary but would come in handy if values were
        # precalculated) we cannot simply continue or there would be
        # some missing. 'continue/next' here if needed.
        if map_to_seq and map_to_seq[cidx] in GAP_CHARS:
            LOG.debug("Skipping col %d because map_to_seq has gap there." % (cidx+1))
            continue

        counts_str = ' '.join(
            ["%s:%d" % (k,v) for (k,v) in sorted(counter.iteritems())])
        if not map_to_seq:
            rep_col = cidx
        else: 
            rep_col = map_to_seq_cols[cidx]
        print "%d %.6f %.6f %s" % (
            rep_col+1 if not map_to_seq else map_to_seq_cols[cidx]+1, 
            seqid_per_col[cidx], entropy_per_col[cidx], counts_str)

    if fh != sys.stdout:
        fh.close()

        
if __name__ == "__main__":
    main()


