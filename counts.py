#!/usr/bin/env python
"""Create histogram of data
"""

#--- standard library imports
#
import sys
import argparse


#--- third-party imports
#
import numpy

#--- project specific imports
#
# /


__author__ = "Andreas Wilm"
__version__ = "0.1"
__email__ = "andreas.wilm@gmail.com"
__license__ = "The MIT License (MIT)"



def main():
    """main function
    """
    
    parser = argparse.ArgumentParser()
    default = 10
    parser.add_argument("--hist",
                        action="store_true",
                        help="Report histogram instead of counts/frequencies")
    parser.add_argument("--bins",
                        type=int,
                        default=default,
                        help="Number of equal-width bins (--hist only)"
                        " (default %d)" % default)
    parser.add_argument("--min",
                        type=float,
                        help="Lower range for bins")
    parser.add_argument("--max",
                        type=float,
                        help="Upper range for bins")
    parser.add_argument("-i", "--infile",
                        default="-",
                        help="Input file containing one value per line"
                        " (- for stdin, which is default)")
    args = parser.parse_args()
    
        
    if args.infile == "-":
        fh = sys.stdin
    else:
        fh = open(args.infile, 'r')
                
    # taken from describe.py: using an iterable seems to be the most
    # efficient way to dynamically grow an array
    #
    iterable = (float(line) for line in fh if len(line.strip())>0)
    arr = numpy.fromiter(iterable, numpy.float)

    if fh != sys.stdin:
        fh.close()


    if args.hist:
        bin_range = (args.min if args.min!=None else arr.min(), 
                     args.max if args.max!=None else arr.max())
        (hist, bin_edges) = numpy.histogram(arr, bins=args.bins, range=bin_range)
        print "#lower bound\tupper bound\tcounts"
        print "#note: all but last (righthand-most) bin are half-open, i.e. [...)"
        for (i, val) in enumerate(hist):
            print "{}\t{}\t{}".format(bin_edges[i], bin_edges[i+1], val)
    else:
        unique, counts = numpy.unique(arr, return_counts=True)
        print "#count\tvalue"
        for (u, c) in zip(unique, counts):
            print "{}\t{}".format(c, u)
        #print numpy.asarray((unique, counts)).T


    

if __name__ == "__main__":
    main()
    
