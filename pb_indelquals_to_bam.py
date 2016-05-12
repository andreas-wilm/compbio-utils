#!/usr/bin/env python
"""Reads indel qualities from binary Pacbio data (BAS/BAX) and adds
them to given BAM (BI/BD GATK-style)
"""


import sys
import os
import logging
import argparse
import pickle
import tempfile

#--- third-party imports
#

    
# global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')


assert (3, 0) > sys.version_info
# python3 not tested and pbcore not installed there anyway

import pysam
pysam_version = tuple([int(x) for x in  pysam.__version__.split('.')])
assert (0, 8, 0) <= pysam_version, (pysam_version)


__author__ = "Andreas Wilm"
__version__ = "0.1"
__email__ = "wilma@gis.a-star.edu.sg"
__license__ = "The MIT License (MIT)"



def phred_to_ascii(p):
    """converted phred score to ascii value"""
    return chr(p+33)


def rotate(l, n):
    """rotate list"""
    return l[-n:] + l[:-n]


def parse_indelquals_from_pb(bas_or_bax, beep_every_x_reads=100):
    """parses indel qualities from Pacbio files and returns them as
    dictionary (read name is key)
    """
    # http://pacificbiosciences.github.io/pbcore/pbcore.io.html
    from pbcore.io import BasH5Reader
    assert (0, 9, 4) == [int(x) for x in pbcore.__VERSION__.split('.')]
    # other versions might work. not tested
    # can't load this and pysam at the same time on our setup because of version conflict

    # WARNING indel quals are read into memory! key = sq
    indelquals = {}

    # BasH5Reader can read bax and bas
    bas = BasH5Reader(bas_or_bax)
    readcount = zmwcount = bascounts = 0
    for bax in bas.parts:
        LOG.info("Parsing %s" % bax.filename)
        bascounts += 1
        for zmw in bax:
            zmwcount += 1
            for r in zmw.subreads:
                readcount += 1
                ins_phred = [phred_to_ascii(x) for x in r.InsertionQV()]
                del_phred = [phred_to_ascii(x) for x in r.DeletionQV()]
                # shift by one to get GATK style: Q that next base is indel error
                ins_phred = ''.join(rotate(ins_phred, -1))
                del_phred = ''.join(rotate(del_phred, -1))
                indelquals[r.readName] = {'BI': "BI:Z:{}".format(ins_phred),
                                          'BD': "BD:Z:{}".format(del_phred)}
                # NOTE: B[ID]:Z: also stored
                if readcount % beep_every_x_reads == 0:
                    LOG.info("{:d} reads parsed...".format(readcount+1))

                #LOG.warn("DEBUG break"); break;
            #LOG.warn("DEBUG break"); break;

    LOG.info("Parsed BI/BD from {:d} reads, {:d} ZMWs and {:d} BAS files".format(
        readcount, zmwcount, bascounts))
    return indelquals


def get_indelquals(bas_or_bax, indelqual_pickle):
    """Extract or load pre-extracted indel qualities"""

    if indelqual_pickle and os.path.exists(indelqual_pickle):
        LOG.info("Loading pre-extracted indel qualities from {}".format(
            indelqual_pickle))
        with open(indelqual_pickle, 'rb') as fh:
            indelquals = pickle.load(fh)
    else:
        if not indelqual_pickle:
            (ofh, indelqual_pickle) = tempfile.mkstemp()
            os.close(ofh)

        indelquals = parse_indelquals_from_pb(bas_or_bax)
        with open(indelqual_pickle, 'wb') as fh:
            pickle.dump(indelquals, fh)
            LOG.info("indel qualities pickled to {}".format(indelqual_pickle))

    return indelquals


def main():
    """main function"""

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("--pacbio",
                        required="True",
                        dest="bas_or_bax",
                        help="Pacbio BAS or BAX to get indel qualities from")
    parser.add_argument("--bam-in",
                        dest="bam_in",
                        required=True,
                        help="BAM to add BI/BD to")
    parser.add_argument("--bam-out",
                        dest="bam_out",
                        required=True,
                        help="BAM to output")
    parser.add_argument("--indelquals",
                        dest="indelqual_pickle",
                        help="Pickle containing already extracted BI/BD values."
                        " If non-existant, write here for later use, otherwise read from here")
    parser.add_argument("--verbose",
                        action="store_true",
                        help="Be verbose")

    args = parser.parse_args()

    if args.verbose:
        LOG.setLevel(logging.INFO)

    for f in [args.bas_or_bax, args.bam_in]:
        assert os.path.exists(f)
    for f in [args.bam_out]:
        assert not os.path.exists(f)
    
    indelquals = get_indelquals(args.bas_or_bax, args.indelqual_pickle)

    
    bam_in = pysam.AlignmentFile(args.bam_in, 'rb')
    bam_out = pysam.AlignmentFile(args.bam_out, "wb", template=bam_in)
    readcount = changedcount = 0
    indelquals_querynames = set(indelquals.keys())# for faster lookup
    for r in bam_in:
        if r.query_name in indelquals_querynames:
            for tag in ['BI', 'BD']:
                q = indelquals[r.query_name][tag]
                if q.startswith("{}:Z:".format(tag)):
                    q = q[5:]
                assert len(q) == len(r.query_sequence)
                r.set_tag(tag, q)
        bam_out.write(r)
        readcount += 1
        changedcount += 1
        if readcount % 100 == 0:
            LOG.debug("Written {} reads ({} with updated BI/BD)".format(
                readcount, changedcount))

if __name__ == "__main__":
    main()
