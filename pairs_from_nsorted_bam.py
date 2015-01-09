#!/usr/bin/env python
"""Extract only paired reads from name sorted BAM files"""

__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2013 Genome Institute of Singapore"
__license__ = "WTFPL"


#--- standard library imports
#
import argparse
import logging
import os
import sys

#--- third-party imports
#
import pysam

#--- project specific imports
#


#--- globals
#
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN, 
                    format='%(levelname)s [%(asctime)s]: %(message)s')



def cmdline_parser():
    """
    creates an argparse instance
    """
    
    # http://docs.python.org/library/optparse.html
    parser = argparse.ArgumentParser(
        description=__doc__)
    
    parser.add_argument("--verbose",
                        dest="verbose",
                        action="store_true",
                        help="Enable verbose output")
    parser.add_argument("--debug",
                        dest="debug",
                        action="store_true", 
                        help=argparse.SUPPRESS) #"debugging")
    parser.add_argument("-i", "--in",
                        dest="bam_in",
                        required=True,
                        help="Input BAM file (must be name sorted)."
                        " '-' for stdin")
    parser.add_argument("-o", "--out",
                        dest="bam_out",
                        required=True,
                        help="Output BAM file. '-' for stdout")
    return parser



def qname_base(qname):
    """Return read name base, i.e. read name with pair information
    """
    
    if qname[-2] in [".", "/", "#"]:
        return qname[:-2]
    else:                                    
        return qname

    
def main():
    """main function
    """

    n_reads_in = 0
    n_reads_out = 0
    
    parser = cmdline_parser()
    args = parser.parse_args()
    
    if args.verbose:
        LOG.setLevel(logging.INFO)
    if args.debug:
        LOG.setLevel(logging.DEBUG)
            
    if not os.path.exists(args.bam_in) and args.bam_in != "-":
        LOG.error("File '%s' does not exist.\n" % (args.bam_in))
        sys.exit(1)                                            
    if os.path.exists(args.bam_out) and args.bam_out != "-":
        LOG.error("Cowardly refusing to overwrite file '%s'\n" % (args.bam_out))
        sys.exit(1)                                            
        
    sam_in = pysam.Samfile(args.bam_in, "rb" )
    sam_out = pysam.Samfile(args.bam_out, "wb", template=sam_in)

    last_read = None
    for read in sam_in:
        n_reads_in += 1
        
        if last_read == None:
            #print "First %s" % read.qname
            last_read = read
        elif qname_base(last_read.qname) == qname_base(read.qname):
            #print "Writing pair %s %s" % (last_read.qname, read.qname);
            sam_out.write(last_read)
            sam_out.write(read)
            n_reads_out += 2
            last_read = None
        else:
            #print "Dropping %s" % last_read.qname
            #print "First %s" % read.qname
            last_read = read

    LOG.info("%d reads in. %d reads out." % (n_reads_in, n_reads_out))
    sam_in.close()
    sam_out.close()


    
if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
        
