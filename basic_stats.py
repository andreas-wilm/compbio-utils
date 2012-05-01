#!/usr/bin/env python
"""Compute some basic statistics on numerical data, read as one data
point per line.
"""

#--- standard library imports
#
import sys

#--- third-party imports
#
import numpy

#--- project specific imports
#
# /


__author__ = "Andreas Wilm"
__version__ = "0.1"
__email__ = "andreas.wilm@gmail.com"
__copyright__ = ""
__license__ = ""
__credits__ = [""]
__status__ = ""


def basic_stats(arr):
    """Compute basic statistics on the given data
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
        
    

if __name__ == "__main__":
    
    try:
        fin = sys.argv[1]
        if fin == "-":
            fh = sys.stdin
        else:
            fh = open(fin, 'r')
    except IndexError:
        sys.stderr.write("Expecting file name as input. Use '-' for stdin\n")
        sys.exit(1)

    # using an iterable seems to be the most efficient way to
    # dynamically grow an array
    iterable = (float(line) for line in fh)
    arr = numpy.fromiter(iterable, numpy.float)

    res_list = basic_stats(arr)
    
    for (name, res) in res_list:
        print "%s\t%f" % (name, res)

    if fh != sys.stdin:
        fh.close()
