#!/usr/bin/env python

import glob
import gzip
import os, sys
from time import strftime
#import itertools

import numpy
#from scipy.sparse import lil_matrix

def now():
    return strftime("%Y-%m-%d %H:%M:%S")


def parse_rnaplfold(bp_or_dp_file):
    """parse rnaplfold plain basepair or dp.ps file (gzip supported)

    format of dp.ps
    ---------------
    %This file contains the square roots of the base pair probabilities in the
    % form
    % i  j  sqrt(p(i,j)) ubox

    basepair files are just
    -----------------------
    i  j  sqrt(p(i,j))

    Returns symmetric, 2D-array as dict where keys are zero-offset
    positions i and j, with i<j and max pos seen (also zero-offset)
    """

    bp_mat = dict()

    if bp_or_dp_file[-3:] == '.gz':
        fh = gzip.open(bp_or_dp_file)
    else:
        fh = open(bp_or_dp_file)

    max_pos = -1
    for line in fh:
        line_split = line.rstrip().split()

        # dp.ps
        if len(line_split)==4 and not line_split[3]=="ubox":
            #print "DEBUG: ignoring:", line
            continue
        elif len(line_split)==3 and not line_split[0].isdigit():
            #print "DEBUG: ignoring:", line
            continue

        i = int(line_split[0])-1
        j = int(line_split[1])-1
        sqrt_p = float(line_split[2])

        assert i < j
        key = (i, j)
        assert bp_mat.has_key(key) == False
        bp_mat[key] = sqrt_p**2

        max_pos = max(j, max_pos)
        # max j might not be at last i
    fh.close()

    return (bp_mat, max_pos)


def main():
    """main function
    """

    bp_files = glob.glob('genome_rnaplfold_dp_ps/*basepairs.gz')
    #sys.stderr.write("FIXME: fixing file to get small one\n"); bp_files = ['genome_rnaplfold_dp_ps/chrI_basepairs.gz']

    wig_file_up = "rnaplfold_vector_up.wig"
    wig_file_down = "rnaplfold_vector_down.wig"
    wig_file_merged = "rnaplfold_vector_merged.wig"
    for f in [wig_file_up, wig_file_down, wig_file_merged]:
        assert not os.path.exists(f), ("Cowardly refusing to overwrite %s" % f)
    fh_out_up = open(wig_file_up, 'w')
    fh_out_down = open(wig_file_down, 'w')
    fh_out_merged = open(wig_file_merged, 'w')
    
    for f in bp_files:
        chrom = os.path.basename(f).replace("_basepairs.gz", "")
        sys.stderr.write("INFO(%s): parsing %s\n" % (now(), f))
        (bp_dict, max_pos) = parse_rnaplfold(f)
        sys.stderr.write("INFO(%s): got %d values with max_pos %d\n" % (now(), len(bp_dict), max_pos))
        # create a max_pos by max_diff array instead of max_pos to
        # max_pos to save memory (and time)
        #max_diff = max([x[1]-x[0] for x in bp_dict.keys()])
        #bp_mat = numpy.zeros((max_pos, max_diff))

        #print "INFO(%s): declaring matrix" % (now())
        #bp_mat = lil_matrix((max_pos+1, max_pos+1))
        # how to efficiently iterate over a sparse matrix
        # http://stackoverflow.com/questions/4319014/iterating-through-a-scipy-sparse-vector-or-matrix
        #for (i, j, v) in itertools.izip(bp_mat.row, bp_mat.col, bp_mat.data):

        # up is paired upstream i.e. towards 5'. down is reverse
        up_arr = numpy.zeros(max_pos+1)
        down_arr = numpy.zeros(max_pos+1)
        mod = len(bp_dict.items())/20
        for (i, (k, v)) in enumerate(bp_dict.items()):
            if i%mod==0:
                sys.stderr.write("INFO(%s): done with converting %d%% of values to proper matrix\n" % (
                    now(), (i+1)/float(len(bp_dict))*100.0))
            #bp_mat[k[0], k[1]] = v
            #bp_mat[k[1], k[0]] = v
            assert k[0] < k[1]
            down_arr[k[0]] += v
            up_arr[k[1]] += v
        #import pdb; pdb.set_trace()
        
        #http://genome.ucsc.edu/goldenPath/help/wiggle.html
        merged_arr = up_arr + down_arr
        for (ori, arr, fh_out) in [('up', up_arr, fh_out_up),
                                   ('down', down_arr, fh_out_down),
                                   ('overall', merged_arr, fh_out_merged)]:
            fh_out.write('track type=wiggle_0 name="bpp-%s-%s" description="RNAplfold %s Pairing Vector %s" visibility=full autoScale=on\n' % (ori, chrom, ori, chrom))
            fh_out.write("variableStep chrom=%s\n" % (chrom))
            it = numpy.nditer(arr, flags=['c_index'])
            while not it.finished:
                if it[0]>0.000000:
                    fh_out.write("%d %f\n" % (it.index+1, it[0]))
                it.iternext()
    for fh_out in [fh_out_up, fh_out_down, fh_out_merged]:
        fh_out.close()

if __name__ == "__main__":
    main()
