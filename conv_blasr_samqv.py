#!/usr/bin/env python
"""Converts blasr -printSAMQV SAM format as follows:
- use subst qual for base qual
- use qi as BI with offset
- use qd as BD with offset
"""

__author__ = "Andreas Wilm"
__version__ = "0.1"
__email__ = "wilma@gis.a-star.edu.sg"
__license__ = "The MIT License (MIT)"


import sys
import os
from collections import OrderedDict

import pysam


def conv_read(r):
    # in-place editing of tags and qualities apparently not possible in 
    # pysam (see also FAQ)
        
    tags = r.tags
    # using a dict for convenience. also gets rid of ts duplicate
    dtags = OrderedDict(tags)
    # offset needed to make it BI and BD which are quals after current 
    # base: take 1+ and add dummy at end
    bi = dtags['qi'][1:] + '#'
    bd = dtags['qd'][1:] + '#'
    dtags['BI'] = bi
    dtags['BD'] = bd


    # pacbio's base qualities are mergers of all qualities. 
    # we want subst quals instead
    if True:
        dtags['qo'] = r.qual# save original
        qual = r.qual
        qual = dtags['qs']
        r.qual = qual
    
    r.tags = dtags.items()
    # FIXME pysam bug: rg:z: to rg:a:
    # use set_tag instead

def main(sam_in, sam_out):
    """main function"""
    
    for r in sam_in:
        
        #print "BEFORE"
        #sam_out.write(r)
        conv_read(r)
        #print "AFTER"
        sam_out.write(r)
        #sys.stderr.write("DEBUG exit\n"); sys.exit(1)
        
    
    
if __name__ == "__main__":
    assert len(sys.argv)==3, ("Usage: %s basr-samqv-in.bam conv-out.bam" % (
        os.path.basename(sys.argv[0])))

    samfile_in = sys.argv[1]
    sam_in = pysam.Samfile(samfile_in)
    
    samfile_out = sys.argv[2]
    if samfile_out == "-":
        mode = "wb"
    else:
        assert not os.path.exists(samfile_out)
        mode = "wb"
    sam_out = pysam.Samfile(samfile_out, mode, template=sam_in)

    main(sam_in, sam_out)
    
    sam_in.close()
    sam_out.close()
