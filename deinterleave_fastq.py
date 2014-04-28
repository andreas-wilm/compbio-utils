#!/usr/bin/env python
"""Deainterleave FastQ. Can handle paired and single end.
"""

__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2014 Genome Institute of Singapore"
__license__ = "GPL2"


#--- standard library imports
#
import sys
import logging
import os
import gzip
from collections import namedtuple

#--- third-party imports
#
#/

#--- project specific imports
#
# /


FastqSeq = namedtuple('FastqSeq', ['id', 'sq', 'bq'])


#global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')


def get_baseid(seqid):
    """return base of sequence id, i.e. without pair info.

    would also work if we used biopython ids, because for newer
    fastq (illumina>1.8+) the pair info is in the description not the
    id/name so there it's safe to return id as it is
    """

    if seqid[-2] == "/":
        #print "returning %s for %s" %  (seqid[:-2], seqid)
        return seqid[:-2]
    else:
        #print "returning %s" %  seqid
        return seqid


def parse_fastq(fh):
    """iterator for reading fastq file. keep it simple and make no
    interpretation of data, e.g. base qualities (exception is id2
    which is not saved). this way we save some time and each entry can
    almost be printed as is"""

    while True:
        # check for eof which should only happen after complete entry,
        # i.e before new id. "readline(): empty string is returned
        # only when EOF is encountered immediately."
        sid = fh.readline()
        if len(sid)==0:
            break
        assert sid[0] == '@', (
            "Expected id starting with @, but got '%s'" % sid)
        sid = sid.rstrip()[1:]

        sq = fh.readline().rstrip()
        tmp = fh.readline().rstrip()
        assert tmp[0] == "+", (
            "Expected second id line starting with +, but got '%s'" % tmp)
        bq = fh.readline().rstrip()
        assert len(sq) == len(bq), (
            "Length mismatch between bases (%d) and base qualities (%d)" % (len(sq), len(bq)))

        yield FastqSeq(sid, sq, bq)


def write_fastq(fh, fastqseq):
    fh.write("@%s\n%s\n+\n%s\n" % (fastqseq.id, fastqseq.sq, fastqseq.bq))


def main():
    """main function
    """

    # poor man's argparse
    try:
        fq_in = sys.argv[1]
        fq_outbase = sys.argv[2]
    except IndexError:
        sys.stderr.write("Usage: %s in.fq[.gz] out_fq_prefix\n" % sys.argv[0])
        sys.exit(1)
        
    LOG.setLevel(logging.INFO)

    if fq_in == "-":
        fh_in = sys.stdin
    else:
        assert os.path.exists(fq_in)
        if fq_in[-3:] == ".gz":
            fh_in = gzip.open(fq_in)
        else:
            fh_in = open(fq_in)

    fq_out_pe1 = fq_outbase + "_1.fastq.gz"
    fq_out_pe2 = fq_outbase + "_2.fastq.gz"
    fq_out_sr = fq_outbase + "_0.fastq.gz"
    for f in [fq_out_pe1, fq_out_pe2, fq_out_sr]:
        if os.path.exists(f):
            LOG.fatal("Cowardly refusing to overwrite already existing file %s" % f)
            sys.exit(1)
    LOG.info("Writing PE reads to %s and %s and single end reads to %s" % (
        fq_out_pe1, fq_out_pe2, fq_out_sr))

    fh_out_pe1 = gzip.open(fq_out_pe1, 'w')
    fh_out_pe2 = gzip.open(fq_out_pe2, 'w')
    fh_out_sr = gzip.open(fq_out_sr, 'w')

    prevseq = None
    curseq = None
    num_in = num_sr = num_pe = 0
    for curseq in parse_fastq(fh_in):
        num_in += 1
        if prevseq:
            # name match indicates pair, otherwise the previously read seq was single
            if get_baseid(prevseq.id)==get_baseid(curseq.id):
                write_fastq(fh_out_pe1, prevseq)
                write_fastq(fh_out_pe2, curseq)
                prevseq = None
                num_pe += 2
            else:
                write_fastq(fh_out_sr, prevseq)
                num_sr += 1
                prevseq = curseq
        else:
            prevseq = curseq
    # if there's still on left, then it must have been sr
    if prevseq:
        write_fastq(fh_out_sr, prevseq)
        num_sr += 1

    for fh in [fh_out_pe1, fh_out_pe2, fh_out_sr, fh_in]:
        fh.close()

    LOG.info('Written %d PE reads and %d SR reads (total was %d)' % (
        num_pe, num_sr, num_in))


if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
