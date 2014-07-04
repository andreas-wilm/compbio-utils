#!/usr/bin/env python
"""Performs a mapping of given SR/PE reads (gzip supported) with
BWA-MEM against given source of contamination and produces an
(unsorted) BAM file with contaminated reads (one mate mapping suffices
to make pair count as contamination) and new, gzipped fastq file(s)
with clean reads (qualities will be Sanger encoded!).

Needs samtools and bwa (mem) installed.
"""

__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2014 Genome Institute of Singapore"
__license__ = "WTFPL http://www.wtfpl.net/"


#--- standard library imports
#
import sys
import logging
import os
import argparse
import subprocess
from collections import deque
from string import maketrans
import gzip

#--- third-party imports
#
#import pysam


# global logger
#
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')


def read_base_name(r):
    if r[-1] in "0123456789" and r[-2] in "#/":
        return r[:-2]
    else:
        return r


# http://stackoverflow.com/questions/1738633/more-pythonic-way-to-find-a-complementary-dna-strand
def complement(strand):
    return strand.translate(maketrans('TAGCtagc', 'ATCGATCG'))


def sam_to_fastq(sam_line, fastq_fh, check_uniq_occurance=10000):
    """convert sam line to fastq entry.

    Will make an attempt to check that no reads with identical names
    were processed (keeping check_uniq_occurance entries). The larger
    the number the more things will slow down. This makes it in fact
    the most time consuming effect, which should be unnecessary as
    long as we pass 0x900==0 reads in here.
    """

    # local static fake
    # http://stackoverflow.com/questions/279561/what-is-the-python-equivalent-of-static-variables-inside-a-function
    if 'seen' not in sam_to_fastq.__dict__:
        sam_to_fastq.seen = deque(maxlen=check_uniq_occurance)

    line_split = sam_line.split('\t')
    (name, flag) = (line_split[0:2])
    name = read_base_name(name)# necessary?
    flag = int(flag)
    (seq, qual) = line_split[9:11]

    if flag & 0x10:
        seq = reversed(complement(seq))
        qual = reversed(qual)

    # From SAM spec: "If 0x1 is unset, no assumptions can be made
    # about 0x2, 0x8, 0x20, 0x40 and 0x80.". So check for pair before
    # checking read number.
    index = 1
    if (flag & 0x1) and (flag & 0x80):
        index += 1
    name = "%s/%d" % (name, index)

    if check_uniq_occurance:
        assert name not in sam_to_fastq.seen
        sam_to_fastq.seen.append(name)

    fastq_fh.write('@%s\n%s\n+\n%s\n' % (name, seq, qual))


def main(fastq_in, ref, fastq_fh, bam_out, num_threads=2):
    """main function

    fastq_in and fastq_fh should be lists (single or paired-end)
    """

    # bufsize needs testing & optimization
    bufsize = 4096

    assert len(fastq_in) == len(fastq_fh)
    for f in fastq_in + [ref]:
        assert os.path.exists(f)

    if len(fastq_in) == 2:
        last_mate_name = None

    bwa_mem_cmd = ['bwa', 'mem', '-t', "%d" % num_threads, ref]
    bwa_mem_cmd.extend(fastq_in)

    # can't use pysam with subprocess
    #p = subprocess.Popen(bwa_mem_cmd, stdout=subprocess.PIPE)
    #s = pysam.Samfile(p.stdout, "rb")
    # throws 'TypeError: Argument must be string or unicode'
    # see also https://www.biostars.org/p/15298/#105456
    # could use stringio?

    cmd = ["samtools", "view", "-bS", "-"]
    samtools_p = subprocess.Popen(
        cmd, stdin=subprocess.PIPE, stdout=open(bam_out, 'wb'), bufsize=bufsize)

    bwa_p = subprocess.Popen(bwa_mem_cmd,
                             stdout=subprocess.PIPE, bufsize=bufsize)
    for line in bwa_p.stdout:
    #for (line_no, line) in enumerate(bwa_p.stdout):
        #if line_no>1000:
        #    LOG.critical("DEBUG break")
        #    break

        # catch sam header
        if line.startswith('@'):
            samtools_p.stdin.write(line)
            continue

        (name, flag) = line.split('\t')[:2]
        flag = int(flag)
        name = read_base_name(name)# necessary?

        # from the SAM spec:
        # http://samtools.github.io/hts-specs/SAMv1.pdf:
        #
        # - Bit 0x4 is the only reliable place to tell whether the
        # read is unmapped. If 0x4 is set, no assumptions can be made
        # about RNAME, POS, CIGAR, MAPQ,
        #
        # - If 0x1 is unset, no assumptions can be made about 0x2,
        # 0x8, 0x20, 0x40 and 0x80.
        #
        # - For each read/contig in a SAM file, it is required that
        # one and only one line associated with the read satisfies
        # 'FLAG & 0x900 == 0'. This line is called the primary line of
        # the read.
        #
        # - Bit 0x100 marks the alignment not to be used in certain
        # analyses when the tools in use are aware of this bit. It is
        # typically used to flag alternative mappings when multiple
        # mappings are presented in a SAM.
        #
        # - Bit 0x800 indicates that the corresponding alignment line
        # is part of a chimeric alignment. A line flagged with 0x800 is
        # called as a supplementary line

        if flag & 0x900:
            LOG.debug("Skipping secondary or"
                      " supplementary alignment: %s" % line)
            continue

        is_paired = flag & 0x1
        if len(args.fastq_in) == 2:
            assert is_paired, (
                "Got two fastq files (paired end) but bwa mem output"
                 " reported single end mapping")
        else:
            assert not is_paired, (
                "Got one fastq files (single end) but bwa mem output"
                 " reported paired end mapping")

        if is_paired:
            if last_mate_name:
                assert name == last_mate_name
                last_mate_name = None
            else:
                last_mate_name = name

            both_unmapped = (flag & 0x4) and (flag & 0x8)
            if not both_unmapped:
                samtools_p.stdin.write(line)
            else:
                if flag & 0x40:
                    sam_to_fastq(line, fastq_fh[0])
                elif flag & 0x80:
                    sam_to_fastq(line, fastq_fh[1])
                else:
                    raise ValueError()
        else:
            unmapped = flag & 0x4
            if unmapped:
                sam_to_fastq(line, fastq_fh[0])
            else:
                samtools_p.stdin.write(line)


    samtools_p.stdin.close()
    if samtools_p.wait() != 0:
        LOG.critical("unhandled samtools error (this can happen"
                     " if no read counted as contaminated)."
                     " You better tally numbers manually")



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # FIXME allow for bwa-mem special args
    # FIXME allow for reuse of existing bam (needs to be sorted by name!)
    #parser.add_argument("-b", "--bam",
    #                    dest='bam',
    #                    help="Already existing BAM file mapped against source of contamination")
    mandatory = parser.add_argument_group('mandatory arguments')
    mandatory.add_argument("-i", "--fq",
                        dest='fastq_in',
                        required=True,
                        nargs="*",
                        help="FastQ input file/s")
    mandatory.add_argument("-o", "--outpref",
                        dest='outpref',
                        required=True,
                        help="Filename prefix for output files")
    mandatory.add_argument("-r", "--ref",
                        required=True,
                        dest='ref',
                        help="Reference fasta file of source of contamination"
                        " (needs to be bwa indexed already)")
    default = 2
    mandatory.add_argument("-t", "--threads",
                        dest='num_threads',
                        default=default,
                        help="Number of threads to use for mapping"
                        " (default = %d)" % default)
    args = parser.parse_args()

    fastq_out = ["%s_%d.fastq.gz" % (args.outpref, i+1)
                 for i in range(len(args.fastq_in))]
    bam_out = "%s.bam" % (args.outpref)

    infiles = args.fastq_in + [args.ref]
    for f in infiles:
        if not os.path.exists(f):
            LOG.fatal("Input file %s does not exist" % f)
            sys.exit(1)

    outfiles = fastq_out + [bam_out]
    for f in outfiles:
        if os.path.exists(f):
            LOG.fatal("Cowardly refusing to overwrite"
                      " already existing file %s" % f)
            sys.exit(1)

    fastq_fh = [gzip.open(f, 'w') for f in fastq_out]

    if not os.path.exists(args.ref + ".bwt"):
        LOG.warn("Doesn't look like reference was indexed."
                 " Did you forget to run bwa index %s ?" % args.ref)

    main(args.fastq_in, args.ref, fastq_fh, bam_out,
         args.num_threads)

    for f in fastq_fh:
        f.close()

    LOG.info("Successful program exit")
