#!/usr/bin/env python
"""Decode ASCII-encoded phred scores.

Can replace certain columns with decoded values and leave the rest. Try e.g.
samtools mpileup -s your.bam  | phred-ascii-to-int.py -c 6 7
"""


#--- standard library imports
#
import sys
import argparse

#--- third-party imports
#
# /

#--- project specific imports
#
# /                                                    


__author__ = "Andreas Wilm"
__version__ = "0.1"
__email__ = "andreas.wilm@gmail.com"
__license__ = "The MIT License (MIT)"


def decode_phred(ascii_enc_phred, offset=33):
    """Decode ASCII encoded phred qualities in ascii_enc_phred and
    yield decoded values as int
    """

    for c in ascii_enc_phred:
        res = ord(c)-offset
        assert res >= 0 and res <= 100, ("Phred quality for %c out of range" % c)
        yield res


def main():
    """The main function
    """

    # http://docs.python.org/dev/howto/argparse.html
    parser = argparse.ArgumentParser()
    #parser.add_argument("-v", "--verbose",  action="store_true", # action="count",  default=0,
    #                    help="increase output verbosity")
    default = 33
    parser.add_argument("-e", "--qualenc", 
                        choices=[33, 64], type=int, default=default,
                        help="Qualities are ASCII-encoded Phred +33"
                        " (e.g. Sanger, SRA, Illumina 1.8+) or +64"
                        " (e.g. Illumina 1.3-1.7). Default: %d" % default)
    parser.add_argument("-c", "--col", 
                        type=int, nargs='*',
                        help="Decode qualities in this column"
                        " (tab delimited fields) and replace contents")
    args = parser.parse_args()
    
    # FIXME add file option
    fh = sys.stdin

    for line in fh:
        line = line.rstrip()
        if not len(line):
            print
            continue
        if not args.col:
            enc_quals = line
            dec_quals = decode_phred(enc_quals, args.qualenc)
            print ' '.join(["%s" % q for q in dec_quals])
        else:
            for (idx, col) in enumerate(line.split('\t')):
                if idx > 0:
                    print "\t",
                if idx+1 not in args.col:
                    print col,
                else:
                    enc_quals = col
                    dec_quals = decode_phred(enc_quals, args.qualenc)
                    print ' '.join(["%s" % q for q in dec_quals]),
        print                    
    if fh != sys.stdin:
        fh.close()
                    
                   
if __name__ == "__main__":
    main()
                            
