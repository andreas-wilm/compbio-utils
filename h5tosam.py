#!/usr/bin/env python
"""Converts Pacbio reads to (unaligned) SAM and includes indel qualities (BI 
and BD tag). If BAS files is given all BAX files referenced will be converted.
Otherwise a BAX file is expected.
"""
# TODO
# - multiprocessing option


__author__ = "Andreas Wilm"
__version__ = "0.1"
__email__ = "andreas.wilm@gmail.com"
__license__ = "The MIT License (MIT)"


import sys
import os
from datetime import datetime

# http://pacificbiosciences.github.io/pbcore/pbcore.io.html
from pbcore.io import BasH5Reader


def log(msg):
    """write log message incl timestamp and added newline"""
    t = datetime.isoformat(datetime.now())
    sys.stderr.write("%s %s\n" % (t, msg))

    
def phred_to_ascii(p):
    """converted phred score to ascii value"""
    return chr(p+33)


def rotate(l, n):
    """rotate list"""
    return l[-n:] + l[:-n]

def read_to_sam(r):
    """convert subread to sam entry and return as string"""

    assert len(r.basecalls()) == len(r.SubstitutionQV()) == \
      len(r.DeletionQV()) == len(r.InsertionQV())

    flag = "%d" % 4
    rname = '*'
    pos = "%d" % 0
    mapq = "%d" % 255
    cigar = '*'
    rnext = '*'
    pnext = "%d" % 0
    tlen = "%d" % 0
    log("WARNING: unclear if indel qualities are off by one with respect to GATK (before after current base")
    # From http://files.pacb.com/software/smrtanalysis/2.0.0/doc/smrtportal/help/Portal_PacBio_Glossary.htm:
    # Quality Value (QV): The total probability that the basecall is an insertion or substitution or is preceded by a deletion.
    # Insertion QV: The probability that the basecall is an insertion with respect to the true sequence.
    # Deletion QV: The probability that a deletion error occurred before the current base.
    # Substitution QV: The probability that the basecall is a substitution.
    subst_phred = ''.join([phred_to_ascii(x) for x in r.SubstitutionQV()])
    del_phred = ''.join([phred_to_ascii(x) for x in r.DeletionQV()])
    ins_phred = ''.join([phred_to_ascii(x) for x in r.InsertionQV()])
    bi = "BI:Z:%s" % (rotate(ins_phred, -1))
    bd = "BD:Z:%s" % (rotate(del_phred, -1))
    
    return "\t".join([r.readName, flag, rname, pos, mapq, cigar, 
                     rnext, pnext, tlen, r.basecalls(), subst_phred, bi, bd])


def main(filename, beep_every_x_reads=100):
    """main function"""

    # BasH5Reader can read bax and bas
    bas = BasH5Reader(filename)
    rcount = 0
    print "@HD\tVN:1.5\tSO:unknown"
    for bax in bas.parts:
        log("Parsing %s" % bax.filename)
        for zmw in bax:
            for r in zmw.subreads:
                print read_to_sam(r)
                rcount += 1
                if rcount % beep_every_x_reads == 1:
                    log("%d reads parsed..." % (rcount+1))
    log("Succesful exit after writing %d reads..." % (rcount+1))
                
    
if __name__ == "__main__":
    assert len(sys.argv)==2, ("Usage: %s bas-or-bax" % (
        os.path.basename(sys.argv[0])))
    main(sys.argv[1])
    
