#!/usr/bin/env python
"""Determine number of reads in (gzipped) fastq files
"""


import sys
import os
import gzip
import io


def warn(msg):
    """prints warning message to stderr (appends newline)"""
    sys.stderr.write("WARN: {}\n".format(msg))


def fastqdata_from_fastq(fastq):
    """Check whether we have a data file produced with fastqc for this
    fastq file
    """
    fastqc_dir = fastq
    for ext in [".gz", ".fastq", ".fq"]:
        if fastqc_dir.endswith(ext):
            fastqc_dir = fastqc_dir[:-len(ext)]
    fastqc_dir += "_fastqc"
    fastqc_data = os.path.join(fastqc_dir, "fastqc_data.txt")
    if os.path.exists(fastqc_data):
        if os.path.getctime(fastqc_data) < os.path.getctime(fastq):
            warn("Using fastqc file that's older than fastq file: {}.".format(
                    fastqc_data))
        return fastqc_data
    else:
        return None


def fastq_num_reads(fastq):
    """Determine number of reads in fastq file"""
    fastqc_data = fastqdata_from_fastq(fastq)
    if fastqc_data:
        #sys.stderr.write("DEBUG: is fastqc\n")
        with open(fastqc_data) as fh:
            for line in fh:
                if line.startswith("Total Sequences"):
                    return int(line.split()[-1])
        raise ValueError(fastqc_data)
    else:
        #sys.stderr.write("DEBUG: count fastq.gz\n")
        if fastq.endswith(".gz"):
            # gzip speedup with io.BufferedReader
            # see http://aripollak.com/pythongzipbenchmarks/
            # and http://www.reddit.com/r/Python/comments/2olhrf/fast_gzip_in_python/
            fh = io.BufferedReader(gzip.open(fastq))
        else:
            #sys.stderr.write("DEBUG: count fastq\n")
            fh = open(fastq)
        num_lines = 0
        for line in fh:
            num_lines += 1
            if num_lines % 4 == 1:
                # bytes or string?
                if isinstance(line, bytes):
                    c = line.decode()[0]
                else:
                    c = line[0]
                assert c == '@'
        fh.close()
        assert num_lines % 4 == 0
        return num_lines/4


if __name__ == "__main__":
    for f in sys.argv[1:]:
        if not os.path.exists(f):
            warn("non-existant file {}".format(f))
            continue
        print("{}\t{}".format(f, fastq_num_reads(f)))
        
