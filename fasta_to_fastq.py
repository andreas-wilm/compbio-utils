#!/usr/bin/env python
"""
Convert fasta format to fake fastq format
"""

#--- standard library imports
#
import sys
import os
import gzip
import logging
# optparse deprecated from Python 2.7 on
from optparse import OptionParser


#--- third-party imports
#
from Bio import SeqIO

#--- project specific imports
#
# /



__author__ = "Andreas Wilm"
__version__ = "0.1"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = ""
__license__ = ""
__credits__ = [""]
__status__ = ""



# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')




def cmdline_parser():
    """
    creates an OptionParser instance
    """

    # http://docs.python.org/library/optparse.html
    usage = "%prog: convert fasta output to fake fastq format\n" \
            "usage: %prog [options]"
    parser = OptionParser(usage=usage)

    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="be verbose")
    parser.add_option("", "--debug",
                      action="store_true", dest="debug",
                      help="debugging")
    parser.add_option("-i", "--input",
                      dest="ffasta", # type="string|int|float"
                       help="fasta input file")
    parser.add_option("-p", "--pair",
                      dest="pairno", # type="string|int|float"
                      default="1", choices=["1", "2"],
                      help="mate pair number")
    parser.add_option("-o", "--output",
                      dest="ffastq", # type="string|int|float"
                       help="fastq output file")
    return parser



def main():
    """
    The main function
    """
    
    parser = cmdline_parser()
    (opts, args) = parser.parse_args()

    if len(args):
        parser.error("Unrecognized arguments found: %s." % (
            ' '.join(args)))
        sys.exit(1)
        
    if opts.verbose:
        LOG.setLevel(logging.INFO)

    if opts.debug:
        LOG.setLevel(logging.DEBUG)

    if not opts.ffasta:
        parser.error("fasta input file argument missing.")
        sys.exit(1)
    if not os.path.exists(opts.ffasta):
        LOG.fatal(
            "file '%s' does not exist.\n" % opts.ffasta)
        sys.exit(1)

    if not opts.ffastq:
        parser.error("fastq output file argument missing.")
        sys.exit(1)
    if os.path.exists(opts.ffastq):
        LOG.fatal(
            "Refusing to overwrite existing output file '%s'.\n" % (
                opts.ffastq))
        sys.exit(1)


    if opts.ffasta[-3:] == ".gz":
        fhandle_fa = gzip.open(opts.ffasta, 'r')
    else:
        fhandle_fa = open(opts.ffasta, 'r')
        
    if opts.ffastq[-3:] == ".gz":
        fhandle_fq = gzip.open(opts.fastq, 'w')
    else:
        fhandle_fq = open(opts.ffastq, 'w')


    for seqrec in SeqIO.parse(fhandle_fa, "fasta"):
        default_qual = 'h' # have no scores. use highest score seen in example.
        machine = seqrec.id.split()[0]
        lane = 1
        tile = 1
        xpos = 1
        ypos = 1
        pair = "/%s" % opts.pairno
        seq = str(seqrec.seq).upper()
        qual = len(seq) * default_qual

        fastqid = ':'.join([str(x) for x in
                            [machine, lane, tile, xpos, ypos]])
        fastqid = "%s#%s/%s" % (fastqid, 0, opts.pairno)
        # index is the number for a multiplexed sample (0 for no indexing)
        # index is the barcode in qiime format

        fhandle_fq.write("@%s\n" % fastqid)
        fhandle_fq.write("%s\n" % seq)
        fhandle_fq.write("+%s\n" % fastqid)
        fhandle_fq.write("%s\n" % (default_qual*len(seq)))

    fhandle_fa.close()
    fhandle_fq.close()


     
    
if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
