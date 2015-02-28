#!/usr/bin/env python
"""Split name sorted PE bam randomly into two files keeping pairing intact

create input with like:
$ samtools sort -n -@ 8 coord-sorted.bam name-sorted
"""


#--- standard library imports
#
import sys
import os

# --- third party imports
#
import pysam


#--- project specific imports
#
# /


__author__ = "Andreas Wilm"
__version__ = "0.1"
__email__ = "andreas.wilm@gmail.com"
__license__ = "WTFPL (http://www.wtfpl.net/)"


def main(bam_in, outprefix, splitevery=100000):
    """main function"""

    bam_fh_in = pysam.Samfile(bam_in, "rb")
    bam_fh_out = None
    out_num = -1
    for (n, r) in enumerate(bam_fh_in):
        if n % splitevery == 0:
            if bam_fh_out:
                bam_fh_out.close()
            out_num += 1
            bam_fh_out = pysam.Samfile("%s%d.bam" % (outprefix, out_num),
                                       "wb", template=bam_fh_in)
        bam_fh_out.write(r)
    bam_fh_in.close()

if __name__ == "__main__":
    try:
        bam_in = sys.argv[1]
        outprefix = sys.argv[2]
        splitevery = int(sys.argv[3])
    except IndexError:
        sys.stderr.write("Usage: %s in.bam outprefix splitevery\n", os.path.basename[sys.argv[0]])
        sys.exit(1)
    main(bam_in, outprefix, splitevery)
