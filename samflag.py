#!/usr/bin/env python
"""Explains samflag. Similar to http://picard.sourceforge.net/explain-flags.html
"""

import os, sys

# from http://samtools.sourceforge.net/SAMv1.pdf
# and http://picard.sourceforge.net/explain-flags.html
SAM_FLAG_BITS = [
    (0x1, "template having multiple segments in sequencing", "read paired"),
    (0x2, "each segment properly aligned according to the aligner", "read mapped in proper pair"),
    (0x4, "segment unmapped", "read unmapped"),
    (0x8, "next segment in the template unmapped", "mate unmapped"),
    (0x10, "SEQ being reverse complemented", "read reverse strand"),
    (0x20, "SEQ of the next segment in the template being reversed", "mate reverse strand"),
    (0x40, "the first segment in the template", "first in pair"),
    (0x80, "the last segment in the template", "second in pair"),
    (0x100, "secondary alignment", "not primary alignment"),
    (0x200, "not passing quality controls", "read fails platform/vendor quality checks"),
    (0x400, "PCR or optical duplicate", "read is PCR or optical duplicate"),
    (0x800, "supplementary alignment", "supplementary alignment")]


def explain_samflag(flag):
    """Prints description of active bits in int samflag 
    """

    # If 0x1 is unset, no assumptions can be made about 0x2, 0x8,
    # 0x20, 0x40 and 0x80.
    skip_bits = []
    if not flag & 0x1:
        skip_bits = [0x2, 0x8, 0x20, 0x40, 0x80]
        #print "skip bits %s" % (', '.join(["0x%x" % x for x in skip_bits]))
    for (bit, sam_descr, pic_descr) in SAM_FLAG_BITS:
        if bit in skip_bits:
            continue
        if flag & bit:
            print "0x%x: %s / %s" % (bit, sam_descr, pic_descr)


            
def main():
    """main function
    """
    sys.stderr.write("WARNING: This copy is not maintained anymore. Use the one hosted on https://github.com/CSB5/misc-scripts if in doubt\n")
    
    if '-h' in sys.argv:
        sys.stderr.write("Usage %s [flag]:\n" % (os.path.basename(sys.argv[0])))
        sys.stderr.write("Expects SAM flag as int as only argument."
                         " Otherwise all flags will be listed\n")
        sys.exit(0)
    try:
        flag = int(sys.argv[1])
        explain_samflag(flag)
    except:
        for (bit, sam_descr, pic_descr) in SAM_FLAG_BITS:
            print "0x%x: %s / %s" % (bit, sam_descr, pic_descr)
    

    
        
if __name__ == "__main__":
    main()
