#!/usr/bin/env python
"""Some library routines useful for bioinfo
"""


#--- standard library imports
#
import os

#--- third-party imports
#
#/

#--- project specific imports
#
# /

                                                        
__author__ = "Andreas Wilm"
__version__ = "0.1"
__email__ = "andreas.wilm@gmail.com"
__license__ = "The MIT License (MIT)"


GAP_CHARS = ['-', '~', '.']



def isgap(res):
    """Return true if given residue is a gap character
    """
    return (res in GAP_CHARS)


def ungap(seqstr):
    """Return a copy of sequence string with all gaps removed

    Similar to seq.ungap but aggressively removes all gap characters defined here  
    """
    
    #assert isinstance(seqstr, type("")) or \
    #  isinstance(seqstr, type(u""))
    for c in GAP_CHARS:
        seqstr = seqstr.replace(c, '')
    return seqstr


def guess_seqformat(fseq):
    """Guess sequence format from file extension used by SeqIO
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


