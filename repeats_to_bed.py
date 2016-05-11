#!/usr/bin/env python

# Reads repeat masked (repeats converted to Ns) fasta file and prints
# coordinates of repeat regions in bed format

import os
import sys
from itertools import groupby


def fasta_iter(fasta_name):
    """
    given a fasta file. yield tuples of header, sequence

    Brend Pedersen https://www.biostars.org/p/710/
    """
    fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        #header = header.next()[1:].strip()
        header = next(header)[1:].strip()
        # join all sequence lines to one.
        #seq = "".join(s.strip() for s in faiter.next())
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq

        
def main(fasta, report_repeats=True):
    """main function"""
    for name, seq in fasta_iter(fasta):
        block_start = 0
        for (is_n, grp) in groupby(seq, lambda x: x=='N'):
            l = len(list(grp))
            
            if is_n:
                if report_repeats:
                    print("{}\t{}\t{}".format(name, block_start, block_start+l))
            else:
                if not report_repeats:
                    print("{}\t{}\t{}".format(name, block_start, block_start+l))
            #    print "%s\t%d\t%d" % (seqrec.id, block_start, block_start+l)
            #    # debug
            #    #print seqrec.id, block_start, len(s), ":", s[:5], "...", s[-5:], "==", seqrec.seq[block_start:block_start+5]
            block_start += l

            
def usage():
    """print usage to stderr"""
    myname = os.path.basename(sys.argv[0])
    sys.stderr.write("Usage: {} [-v] seq.fa\n".format(myname))
    sys.stderr.write("{}: report repeat (poly-N) regions in fasta file\n".format(myname))
    sys.stderr.write("Logic is reverted with with -v\n")

                     
if __name__ == "__main__":
    report_repeats = True
    if len(sys.argv) == 2:
        fasta = sys.argv[1]
    elif len(sys.argv) == 3 and sys.argv[1] == "-v":
        report_repeats = False
        fasta = sys.argv[2]
    else:
        usage()
        sys.exit(1)
    assert os.path.exists(fasta)

    main(fasta, report_repeats=report_repeats)
    




















