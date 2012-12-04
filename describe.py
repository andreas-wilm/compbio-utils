#!/usr/bin/env python
"""Compute some basic statistics on numerical data, read as one data
point (float) per line.
"""

#--- standard library imports
#
import sys

#--- third-party imports
#
import numpy
from scipy import stats

#--- project specific imports
#
# /


__author__ = "Andreas Wilm"
__version__ = "0.1"
__email__ = "andreas.wilm@gmail.com"
__license__ = "The MIT License (MIT)"


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
            (numpy.min, None, "minimum"),
            (numpy.max, None, "maximum"),
            (numpy.sum, None, "sum"),
            (numpy.mean, None, "mean"),
            (numpy.std, None, "stdv"),
            (numpy.median, None, "median")
            ]
    for p in [1, 5, 10, 25, 50, 75, 90, 95, 99]:
        func_list.append((numpy.percentile, p, "%dth percentile" % p))

    # call each function and store its name as well as the result
    for (func, arg, name) in func_list:
        if arg:
            result_list.append((name, func(arr, arg)))
        else:
            result_list.append((name, func(arr)))
            
    return result_list
        


def main():
    """
    """
    if len(sys.argv) == 2 and sys.argv[1] != "-":
        fh = open(sys.argv[1], 'r')
    else:
        fh = sys.stdin
                
    # using an iterable seems to be the most efficient way to
    # dynamically grow an array
    #
    iterable = (float(line) for line in fh if len(line.strip())>0)
    arr = numpy.fromiter(iterable, numpy.float)


    # scipy.stats.describe already covers a lot
    describe_res = ["size",
                    "(min, max)", 
                    "arithmetic mean",
                    "unbiased variance",
                    "biased skewness",
                    "biased kurtosis"]
    for (name, res) in zip(describe_res,
                           stats.describe(arr)):
        print "%s:\t%s" % (name, res)
    # ...but not percentiles
    print "median:\t%f" % (stats.scoreatpercentile(arr, 50))
    for p in [1, 5, 10]:
        print "%dth percentile:\t%f" % (p, stats.scoreatpercentile(arr, p))
    upper_q = (stats.scoreatpercentile(arr, 75))
    lower_q = (stats.scoreatpercentile(arr, 25))
    print "IQR:\t%f" % (upper_q-lower_q)

    if fh != sys.stdin:
        fh.close()
    

if __name__ == "__main__":
    main()
    
