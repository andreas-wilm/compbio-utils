#!/usr/bin/env python
"""Unclip soft-clips in BAM file
"""

#--- standard library imports
#
import sys
import os
import argparse
import logging


#--- third-party imports
#
# https://pysam.readthedocs.org/
import pysam
assert [int(x) for x in pysam.__version__.split('.')] >= [0, 8, 2]

# https://github.com/brentp/cigar
import cigar


#--- project specific imports
#
# /


__author__ = "Andreas Wilm"
__version__ = "0.1"
__email__ = "andreas.wilm@gmail.com"
__license__ = "The MIT License (MIT)"


#global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')


def unclip_read(read):
    """unclip read and correct aligned position of necessary"""
    if 'S' in read.cigarstring:
        cig = cigar.Cigar(read.cigarstring)
        offset = 0
        for clen, cop in cig.items():
            if cop != 'S':
                break
            else:
                offset += clen

        cig = cigar.Cigar(read.cigarstring.replace("S", "M"))
        cig = cig.merge_like_ops()
        LOG.debug("Cahnging {}: {} -> {} and pos-= {}".format(read.qname, read.cigarstring, str(cig), offset))
        read.cigarstring = str(cig)
        read.pos -= offset

    return read


def main():
    """main function"""

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--bam-in",
                        required=True, dest='bam_in',
                        help="Input BAM file (to unclip)")
    parser.add_argument("-o", "--bam-out",
                        required=True, dest='bam_out',
                        help="Output BAM file (unclipped)")
    # https://www.reddit.com/r/Python/comments/3nctlm/what_python_tools_should_i_be_using_on_every/
    # by masklinn
    parser.add_argument('-v', '--verbose', action='count', default=0)
    parser.add_argument('-q', '--quiet', action='count', default=0)
    # use as: logging_level = logging.WARN + 10*args.quiet - 10*args.verbose
    # script -vv -> DEBUG
    # script -v -> INFO
    # script -> WARNING
    # script -q -> ERROR
    # script -qq -> CRITICAL
    # script -qqq -> no logging at all

    args = parser.parse_args()

    logging_level = logging.WARN + 10*args.quiet - 10*args.verbose
    LOG.setLevel(logging_level)

    if os.path.exists(args.bam_out):
        LOG.fatal("Cowardly refusing to overwrite existing file {}".format(args.bam_out))
        sys.exit(1)
                 
    LOG.info("Reading from {}".format(args.bam_in))
    LOG.info("Writing to {}".format(args.bam_out))
    # b == binary = bam, i.e. not support for stdin/stdout yet
    samfhin = pysam.AlignmentFile(args.bam_in, 'rb')
    samfhout = pysam.AlignmentFile(args.bam_out, 'wb', template=samfhin)

    num_reads = 0
    for read in samfhin:
        read = unclip_read(read)
        samfhout.write(read)
        num_reads += 1
    LOG.info("{} reads processed".format(num_reads))
    LOG.warn("BAM file might need resorting"
             " and tags incl PNEXT might not be correct anymore\n")

    samfhout.close()
    samfhin.close()


if __name__ == "__main__":
    main()
