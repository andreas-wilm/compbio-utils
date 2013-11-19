#!/usr/bin/env python
"""Compute some basic statistics on numerical data, read as one data
point (float) per line."""

#--- standard library imports
#
import sys

#--- third-party imports
#
import numpy
from scipy.stats import scoreatpercentile, describe

#--- project specific imports
#
# /


__author__ = "Andreas Wilm"
__version__ = "0.1"
__email__ = "andreas.wilm@gmail.com"
__license__ = "The MIT License (MIT)"


def adv_stats(arr):
    """FIXME:add-doc
    """

    print "sum:\t%f" % (arr.sum())
    # scipy.stats.describe already covers a lot
    describe_names = ["size", "(min, max)", "arithmetic mean", 
                      "unbiased variance", "biased skewness", 
                      "biased kurtosis"]
    describe_result = dict(zip(describe_names, describe(arr)))
    for k in describe_names:
        v = describe_result[k]
        if k == "(min, max)":
            # used later
            min = v[0]
            max = v[1]
            continue
        print "%s: %f" % (k, v)
    # ...but not percentiles
    median = scoreatpercentile(arr, 50)
    lower_q = (scoreatpercentile(arr, 25)) # q1
    upper_q = (scoreatpercentile(arr, 75)) # q3
    iqr = upper_q-lower_q
    std = numpy.std(arr)
    print "std:\t%f" % (std)
    #whisker = 1.5*iqr
    #print "IQR:\t%f" % iqr
    for p in [1, 5, 10, 90, 95, 99]:
        print "%dth percentile:\t%f" % (p, scoreatpercentile(arr, p))
    print "# five number summary"
    # FIXME different to http://en.wikipedia.org/wiki/Five-number_summary
    # echo 0, 0, 1, 2, 63, 61, 27, 13, | tr ',' '\n' | describe.py 
    # quartiles differ
    print "min:\t%f" % (min)
    print "q1:\t%f" % (lower_q)
    #print "lower hinge:\t%f" % (median-whisker)
    print "median:\t%f" % (median)
    #print "upper hinge:\t%f" % (median+whisker)
    print "q3:\t%f" % (upper_q)
    print "max:\t%f" % (max)


    
def basic_stats(arr):
    """Compute basic statistics on the given data

    ...
    res_list = basic_stats(arr)
    for (name, res) in res_list:
        print "%s\t%f" % (name, res)

    NOTE: just use scipy.stats.describe instead
    """
    assert arr.ndim == 1
    
    # return
    result_list = []

    # list of functions to call on data, including an extra arg and
    # the functions name
    func_list = [
            (len, None, "length"),
            (numpy.min, None, "min"),
            (numpy.max, None, "max"),
            (numpy.sum, None, "sum"),
            (numpy.mean, None, "mean"),
            (numpy.std, None, "stdv"),
            ]
    for p in [1, 5, 25, 50, 75, 95, 99]:
        if p == 50:
            name = "median"
        elif p == 75:
            name = "q3"
        elif p == 25:
            name = "q1"
        else:
            name = "%dth %%tile" % p
        func_list.append((numpy.percentile, p, name))

    # call each function and store its name as well as the result
    for (func, arg, name) in func_list:
        if arg:
            result_list.append((name, func(arr, arg)))
        else:
            result_list.append((name, func(arr)))
            
    return result_list
        


def main():
    """main function
    """

    VALID_MODES = ['basic', 'advanced']
    
    if len(sys.argv) < 2:
        sys.stderr.write("ERROR: Need at least mode argument (one of %s)\n" % (
            ', '.join(VALID_MODES)))
        sys.exit(1)
        
    mode = sys.argv[1]
    if mode not in VALID_MODES:
        sys.stderr.write("ERROR: First arg needs to be mode (one of %s)\n" % (
            ', '.join(VALID_MODES)))
        sys.exit(1)
        
    if len(sys.argv) == 3 and sys.argv[2] != "-":
        fh = open(sys.argv[2], 'r')
    else:
        fh = sys.stdin
                
    # using an iterable seems to be the most efficient way to
    # dynamically grow an array
    #
    iterable = (float(line) for line in fh if len(line.strip())>0)
    arr = numpy.fromiter(iterable, numpy.float)

    if mode == "basic":
        for (rname, rval) in basic_stats(arr):
            print "%s\t%s" % (rname, rval)
    elif mode == "advanced":
        adv_stats(arr)
    else:
        raise ValueError, ("Unknown mode %s" % mode)
        
    #print "DEBUG: %s" % arr
    if fh != sys.stdin:
        fh.close()
    

if __name__ == "__main__":
    main()
    
