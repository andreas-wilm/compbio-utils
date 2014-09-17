#!/usr/bin/env python
"""Split name sorted PE bam randomly into two files keeping pairing intact

create input with like:
$ samtools sort -n -@ 8 coord-sorted.bam name-sorted
"""


#--- standard library imports
#
import sys
import os
import random
from itertools import izip

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


BEEP_EVERY = 100000

def pairwise(iterable):
    """
    Johnsyweb's answer to Iterating over every two elements in a list
    http://stackoverflow.com/questions/5389507/iterating-over-every-two-elements-in-a-list   
    
    s -> (s0,s1), (s2,s3), (s4, s5), ...
    """
    a = iter(iterable)
    return izip(a, a)


def reads_are_paired(r1, r2, also_check_name=True):
    """ensure two pysam reads are paired
    """
    
    # extra paranoia which might fail if pair no is indicated in name
    if also_check_name and r1.qname != r2.qname:
        return False

    # 0x40: 1st in pair
    # 0x80: 2nd in pair
    # also makes following redundant: assert r1.is_paired and r2.is_paired

    if r1.flag & 0x40:
        return r2.flag & 0x80
    elif r1.flag & 0x80:
        return r2.flag & 0x40
    else:
        return False
            

def main():
    
    f = dict()
    assert len(sys.argv)==5, ("\nUsage: %s in.bam out1.bam out2.bam outs.bam" % (
        os.path.basename(sys.argv[0])))
    f['in'] = sys.argv[1]
    f['out1'] = sys.argv[2]
    f['out2'] = sys.argv[3]
    f['outs'] = sys.argv[4]
    assert os.path.exists(f['in'])

    sam_fh = dict()
    assert os.path.exists(f['in'])
    sam_fh['in'] = pysam.Samfile(f['in'], 'rb')
    num_reads = dict()
    num_reads['in'] = 0
    for k in ['out1', 'out2', 'outs']:
        assert not os.path.exists(f[k])
        sam_fh[k] = pysam.Samfile(f[k], "wb", template=sam_fh['in'])
        num_reads[k] = 0

    for (k, v) in f.items():
        sys.stderr.write("INFO: %s file = %s\n" % (k, v))
        
    for (r1, r2) in pairwise(sam_fh['in']):
        num_reads['in'] += 2
        if num_reads['in'] % BEEP_EVERY == 0:
            for (k, v) in num_reads.items():
                sys.stderr.write("INFO: processed reads for %s: %d\n" % (k, v))

        while not reads_are_paired(r1, r2):
            sys.stderr.write("WARNING: Following two reads don't seem"
                             " to be paired (will write first one to outs"
                             " and try to pair second one): %s and %s\n" % (
                                 r1.qname, r2.qname))
            # write first separately, assuming it's a stray single and next
            # read will belong to second one
            sam_fh['outs'].write(r1)
            num_reads['outs'] += 1
            
            r1 = sam_fh['in'].next()
            num_reads['in'] += 1

        # write
        out_key = random.choice(['out1', 'out2'])
        sam_fh[out_key].write(r1)
        sam_fh[out_key].write(r2)
        num_reads[out_key] += 2

        # debugging only
        #if not random.choice(range(10000)):
        #    sys.stderr.write("DEBUG break leaving truncated output files\n")
        #    break
    
    for k in sam_fh.keys():
        sam_fh[k].close()

    for (k, v) in num_reads.items():
        sys.stderr.write("INFO: processed reads for %s: %d\n" % (k, v))

        
if __name__ == "__main__":
    main()
    sys.stderr.write("INFO: successful exit\n")
    
