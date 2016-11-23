#!/usr/bin/env python

import os
import sys
from collections import namedtuple


Region = namedtuple('Region', ['sq', 'start', 'end'])


MIN_DIST = 1000
MIN_LEN = 10000000


def sq_in_blacklist(sq):
    if sq.startswith("GL"):
        return True
    elif sq.startswith("hs"):
        return True
    elif sq.startswith("NC"):
        return True
    else:
        return False


def region_len(r):
    assert r.start < r.end
    return r.end - r.start


def dist_between_regions(r1, r2):
    """compute distance between two regions on the same
    chromosome. overlap raises error. distance of zero means direct
    neighbours

    FIXME unit tests
    """
    assert r1.sq == r2.sq
    assert r1.start < r1.end
    assert r2.start < r1.end

    if r1.start < r2.start:
        assert r1.end <= r2.start
        d = r2.start - r1.end

    elif r2.start < r1.start:
        assert r2.end <= r1.start
        d = r1.start - r2.end

    else:
        raise ValueError(r1, r2)

    assert d >= 0, (r1, r2)
    return d


def write_region(region, prefix, counter):
    f = "{:s}.{:d}.bed".format(prefix, counter)
    assert not os.path.exists(f)
    with open(f, 'w') as fh:
        fh.write("{:s}\t{:d}\t{:d}\n".format(region.sq, region.start, region.end))


def main(bed, outprefix):
    """main"""

    in_regions = []
    with open(bed) as fh:
        for line in fh:
            ls = line.strip().split()
            sq = ls[0]
            st = int(ls[-2])
            en = int(ls[-1])
            this_region = Region(sq, st, en)
            if sq_in_blacklist(this_region.sq):
                continue
            # check sorting
            if in_regions:
                if in_regions[-1].sq == this_region.sq:
                    assert in_regions[-1].end <= this_region.start
            in_regions.append(this_region)

    last_region = None
    region_out_counter = 0
    for this_region in in_regions:
        if not last_region:
            last_region = this_region

        # new chromosome: print buffer
        elif this_region.sq != last_region.sq:
            write_region(last_region, outprefix, region_out_counter)
            region_out_counter += 1
            last_region = this_region

        # close the previous: join
        elif dist_between_regions(this_region, last_region) < MIN_DIST:
            sys.stderr.write("Region merging. Should have been done with bedtools merge\n")
            last_region = last_region._replace(end=this_region.end)

        # new separate region: old or new one too small
        elif region_len(this_region) < MIN_LEN or region_len(last_region) < MIN_LEN:
            #sys.stderr.write("Merging\n")
            last_region = last_region._replace(end=this_region.end)

        # new separate region and big enough
        else:
            write_region(last_region, outprefix, region_out_counter)
            region_out_counter += 1
            last_region = this_region

    write_region(last_region, outprefix, region_out_counter)
    region_out_counter += 1


if __name__ == "__main__":
    #bed = "gotcloud-bundles_hs37d5-db142-v1_hs37d5.poly-n-split-merge1k.bed"
    #outprefix = "gotcloud-bundles_hs37d5-db142-v1_hs37d5.poly-n-merge1k.split"
    assert len(sys.argv)==3
    bed = sys.argv[1]
    outprefix = sys.argv[2]
    assert os.path.exists(bed)
    
    sys.stderr.write("Reading {}\n".format(bed))
    sys.stderr.write("Writing to {}\n".format(outprefix))
    main(bed, outprefix)
